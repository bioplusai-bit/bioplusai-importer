"""
AlphaMissense import service.
Railway da ishlaydi:
  1. Zenodo dan fayl yuklab oladi
  2. Railway PostgreSQL ga import qiladi
  3. Tugagach o'chadi
"""

import psycopg2
import gzip
import time
import os
import requests
import io

# ── Environment variables (Railway da sozlanadi) ──────────────────────────────
DB_URL    = os.environ.get("DATABASE_URL", "")
ZENODO_URL = "https://zenodo.org/records/10813168/files/AlphaMissense_hg38.tsv.gz"

BATCH_SIZE = 20_000
LOG_EVERY  = 500_000

def main():
    if not DB_URL:
        print("DATABASE_URL environment variable topilmadi!")
        return

    print("=" * 60)
    print("BioPlusAI — AlphaMissense Import Service")
    print("=" * 60)

    print("\n1. PostgreSQL ga ulanmoqda...")
    conn = psycopg2.connect(DB_URL)
    conn.autocommit = False
    cur = conn.cursor()
    print("   OK!")

    # Jadval tekshirish
    cur.execute("""
        SELECT EXISTS (
            SELECT FROM information_schema.tables
            WHERE table_name = 'AlphaMissenseVariants'
        );
    """)
    if not cur.fetchone()[0]:
        print("   AlphaMissenseVariants jadvali topilmadi!")
        print("   Avval migration ishlatib koring.")
        conn.close()
        return

    cur.execute('SELECT COUNT(*) FROM "AlphaMissenseVariants"')
    existing = cur.fetchone()[0]
    print(f"   Hozirgi qatorlar: {existing:,}")

    if existing > 1_000_000:
        print("   Ma'lumotlar allaqachon import qilingan. Chiqilmoqda.")
        conn.close()
        return

    print(f"\n2. Zenodo dan fayl yuklanmoqda...")
    print(f"   URL: {ZENODO_URL}")

    response = requests.get(ZENODO_URL, stream=True, timeout=300)
    response.raise_for_status()

    total_size = int(response.headers.get('content-length', 0))
    print(f"   Hajmi: {total_size / 1024 / 1024:.0f} MB")

    print(f"\n3. Import boshlandi...")
    start_time = time.time()
    total      = 0
    batch      = []
    skipped    = 0
    downloaded = 0

    insert_sql = """
        INSERT INTO "AlphaMissenseVariants"
            ("Chromosome", "Position", "Ref", "Alt", "Score", "Classification", "TranscriptId", "ProteinChange")
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
    """

    # Ustun indekslari
    col_chrom = col_pos = col_ref = col_alt = col_score = col_class = col_trans = col_protein = None
    header_found = False

    # Stream orqali o'qiymiz — to'liq yuklab olmasdan
    buffer = b""
    with requests.get(ZENODO_URL, stream=True, timeout=600) as r:
        r.raise_for_status()

        # gzip stream
        gz = gzip.GzipFile(fileobj=r.raw)

        for raw_line in gz:
            try:
                line = raw_line.decode('utf-8').strip()
            except UnicodeDecodeError:
                continue

            if not line:
                continue

            # Header topish
            if not header_found:
                if 'CHROM' in line.upper() and 'POS' in line.upper():
                    headers = line.lstrip('#').strip().split('\t')
                    headers = [h.strip().lower() for h in headers]
                    print(f"   Header: {headers}")

                    for i, h in enumerate(headers):
                        if h in ('chrom', '#chrom'):      col_chrom   = i
                        elif h == 'pos':                   col_pos     = i
                        elif h == 'ref':                   col_ref     = i
                        elif h == 'alt':                   col_alt     = i
                        elif h == 'am_pathogenicity':      col_score   = i
                        elif h == 'am_class':              col_class   = i
                        elif h == 'transcript_id':         col_trans   = i
                        elif h == 'protein_variant':       col_protein = i

                    header_found = True
                continue

            if line.startswith('#'):
                continue

            parts = line.split('\t')

            try:
                chrom   = parts[col_chrom]
                pos     = int(parts[col_pos])
                ref     = parts[col_ref]   if col_ref   is not None else ''
                alt     = parts[col_alt]   if col_alt   is not None else ''
                score   = float(parts[col_score]) if col_score is not None else 0.0
                cls_raw = parts[col_class] if col_class is not None else ''
                trans   = parts[col_trans]   if col_trans   is not None and len(parts) > col_trans   else None
                protein = parts[col_protein] if col_protein is not None and len(parts) > col_protein else None

                if len(ref) != 1 or len(alt) != 1:
                    skipped += 1
                    continue

                cls = normalize_class(cls_raw)
                batch.append((chrom, pos, ref, alt, score, cls, trans, protein))

                if len(batch) >= BATCH_SIZE:
                    cur.executemany(insert_sql, batch)
                    conn.commit()
                    total += len(batch)
                    batch = []

                    if total % LOG_EVERY == 0:
                        elapsed = time.time() - start_time
                        rate    = total / elapsed
                        eta     = (71_000_000 - total) / rate / 60
                        print(f"   {total:>12,} | {rate:,.0f}/sek | ETA: {eta:.0f} min")

            except (ValueError, IndexError):
                skipped += 1
                continue

    if batch:
        cur.executemany(insert_sql, batch)
        conn.commit()
        total += len(batch)

    elapsed = time.time() - start_time
    print(f"\n4. Import tugadi!")
    print(f"   Jami: {total:,} qator")
    print(f"   Skip: {skipped:,}")
    print(f"   Vaqt: {elapsed:.0f}s ({elapsed/60:.1f} min)")

    # Final count
    cur.execute('SELECT COUNT(*) FROM "AlphaMissenseVariants"')
    final = cur.fetchone()[0]
    print(f"   DB da jami: {final:,} qator")

    conn.close()
    print("\nService tugadi. Railway container o'chadi.")

def normalize_class(raw: str) -> str:
    raw = raw.lower().strip()
    if 'pathogenic' in raw and 'likely' in raw: return 'LIKELY_PATHOGENIC'
    if 'pathogenic' in raw:                     return 'PATHOGENIC'
    if 'benign' in raw and 'likely' in raw:     return 'LIKELY_BENIGN'
    if 'benign' in raw:                         return 'BENIGN'
    if 'ambiguous' in raw:                      return 'AMBIGUOUS'
    return raw.upper()

if __name__ == "__main__":
    main()
