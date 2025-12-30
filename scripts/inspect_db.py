import sqlite3
import pandas as pd
import os


def inspect_database(db_path):
    print(f"\n{'=' * 60}")
    print(f"🕵️‍♂️ 正在体检数据库: {os.path.basename(db_path)}")
    print(f"{'=' * 60}")

    if not os.path.exists(db_path):
        print(f"❌ 文件不存在: {db_path}")
        return

    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # 1. 获取所有表名
        # sqlite_master 是 SQLite 的元数据表，存储了所有对象的信息
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()

        if not tables:
            print("⚠️ 警告: 数据库是空的，没有任何表！")
            conn.close()
            return

        print(f"发现 {len(tables)} 张表:")

        for i, table_tuple in enumerate(tables):
            table_name = table_tuple[0]
            print(f"\n  📄 表 {i + 1}: [{table_name}]")

            # 2. 获取该表的列信息 (PRAGMA table_info)
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = cursor.fetchall()

            # 格式化输出列名和类型
            print(f"     [字段结构]:")
            col_names = []
            for col in columns:
                # col结构: (cid, name, type, notnull, dflt_value, pk)
                c_name = col[1]
                c_type = col[2]
                print(f"       - {c_name} ({c_type})")
                col_names.append(c_name)

            # 3. 检查数据量
            cursor.execute(f"SELECT Count(*) FROM {table_name}")
            row_count = cursor.fetchone()[0]
            print(f"     [数据行数]: {row_count} 行")

            # 4. 预览前 3 行数据
            print(f"     [数据预览]:")
            df_preview = pd.read_sql_query(f"SELECT * FROM {table_name} LIMIT 3", conn)
            print(df_preview.to_string(index=False))
            print("-" * 40)

        conn.close()

    except Exception as e:
        print(f"❌ 无法读取数据库: {e}")


# --- 使用示例 ---
if __name__ == "__main__":
    # 这里填入想检查的文件路径
    inspect_database("./dbs/scrna.db")
    inspect_database("./dbs/atac_breast_cancer.db")
    inspect_database("./dbs/metabolomics.db")
    inspect_database("./dbs/spatial.db")
    inspect_database("./dbs/imaging.db")

