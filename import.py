"""
AlphaMissense import service — MinIO (authenticated) dan o'qiydi.
"""

import psycopg2
import gzip
import time
import os
import io
from minio import Minio

# ── Environment variables ──────────────────────────────────────────────────────
DB_URL         = os.environ.get("DATABASE_URL", "")
MINIO_ENDPOINT = os.environ.get("MINIO_ENDPOINT", "gondola.proxy.rlwy.net:34773")
MINIO_ACCESS   = os.environ.get("MINIO_ACCESS_KEY", "minioadmin")
MINIO_SECRET   = os.environ.get("MINIO_SECRET_KEY", "minioadmin123")
MINIO_BUCKET   = os.environ.get("MINIO_DATASET_BUCKET", "datasets")
MINIO_FILE     = os.environ.get("MINIO_AM_FILE", "AlphaMissense_hg38.tsv.gz")

BATCH_SIZE = 20_000
LOG_EVERY  = 500_000

def main():
    if not DB_URL:
        print("DATABASE_URL topilmadi!")
        return

    print("=" * 60)
    print("BioPlusAI — AlphaMissense Import (MinIO authenticated)")
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
        conn.close()
        return

    cur.execute('SELECT COUNT(*) FROM "AlphaMissenseVariants"')
    existing = cur.fetchone()[0]
    print(f"   Hozirgi qatorlar: {existing:,}")

    if existing > 1_000_000:
        print("   Ma'lumotlar allaqachon import qilingan!")
        conn.close()
        return

    # MinIO ga ulanish
    print(f"\n2. MinIO ga ulanmoqda: {MINIO_ENDPOINT}")
    client = Minio(
        MINIO_ENDPOINT,
        access_key=MINIO_ACCESS,
        secret_key=MINIO_SECRET,
        secure=False
    )

    print(f"   Bucket: {MINIO_BUCKET}")
    print(f"   Fayl: {MINIO_FILE}")

    # Fayl hajmini tekshirish
    try:
        stat = client.stat_object(MINIO_BUCKET, MINIO_FILE)
        print(f"   Hajmi: {stat.size / 1e6:.0f} MB")
    except Exception as e:
        print(f"   Fayl topilmadi: {e}")
        conn.close()
        return

    print(f"\n3. Import boshlandi...")
    start_time   = time.time()
    total        = 0
    batch        = []
    skipped      = 0
    header_found = False

    col_chrom = col_pos = col_ref = col_alt = None
    col_score = col_class = col_trans = col_protein = None

    insert_sql = """
        INSERT INTO "AlphaMissenseVariants"
            ("Chromosome", "Position", "Ref", "Alt", "Score", "Classification", "TranscriptId", "ProteinChange")
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
    """

    try:
        response = client.get_object(MINIO_BUCKET, MINIO_FILE)
        gz = gzip.GzipFile(fileobj=response)

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
                        if h in ('chrom', '#chrom'):  col_chrom   = i
                        elif h == 'pos':               col_pos     = i
                        elif h == 'ref':               col_ref     = i
                        elif h == 'alt':               col_alt     = i
                        elif h == 'am_pathogenicity':  col_score   = i
                        elif h == 'am_class':          col_class   = i
                        elif h == 'transcript_id':     col_trans   = i
                        elif h == 'protein_variant':   col_protein = i

                    header_found = True
                    print(f"   Ustunlar: CHROM={col_chrom}, POS={col_pos}, SCORE={col_score}, CLASS={col_class}")
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

    except Exception as e:
        print(f"   MinIO xato: {e}")
        conn.close()
        return
    finally:
        try:
            response.close()
            response.release_conn()
        except Exception:
            pass

    # Qolgan batch
    if batch:
        cur.executemany(insert_sql, batch)
        conn.commit()
        total += len(batch)

    elapsed = time.time() - start_time
    print(f"\n4. Import tugadi!")
    print(f"   Jami: {total:,} qator")
    print(f"   Skip: {skipped:,}")
    print(f"   Vaqt: {elapsed:.0f}s ({elapsed/60:.1f} min)")

    cur.execute('SELECT COUNT(*) FROM "AlphaMissenseVariants"')
    final = cur.fetchone()[0]
    print(f"   DB da jami: {final:,} qator")

    conn.close()
    print("\nService tugadi!")

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
