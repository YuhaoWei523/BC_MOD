import sqlite3
import pandas as pd

conn = sqlite3.connect("scRNA_data.db")
# 查看前 10 个基因名是什么样子
print("=== 数据库中的前 10 个基因 ===")
df_check = pd.read_sql_query("SELECT DISTINCT Gene FROM Table_Expression LIMIT 10", conn)
print(df_check)

# 模糊搜索 PKM (不管大小写)
print("\n=== 尝试模糊搜索 PKM ===")
df_fuzzy = pd.read_sql_query("SELECT * FROM Table_Expression WHERE Gene LIKE '%PKM%' LIMIT 5", conn)
print(df_fuzzy)

conn.close()