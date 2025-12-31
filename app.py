import asyncio
import sys

# 1. Windows 异步修复
if sys.platform.startswith("win"):
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

import streamlit as st
import sqlite3
import pandas as pd
import os
import json
import zipfile
import io
import time
import base64
import warnings
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from fpdf import FPDF
import auth_manager as auth

# 2. 屏蔽非必要的 UserWarning
warnings.filterwarnings("ignore", category=UserWarning)

# ==========================================
# 3. 国际化字典 (i18n) - 全面覆盖版
# ==========================================
TRANS = {
    "CN": {
        "title": "BC-MOD 乳腺癌多组学数据库",
        "sidebar_header_user": "👤 用户中心",
        "sidebar_header_nav": "🧭 功能导航",
        "sidebar_header_settings": "⚙️ 设置与工具",
        "nav_query": "🔍 数据查询",
        "nav_admin": "🛠️ 后台管理",
        "login_title": "🔐 系统登录",
        "login_btn": "登录",
        "logout_btn": "退出登录",
        "intro_title": "系统介绍",
        "intro_text": "本系统整合 scRNA-seq, ATAC-seq, Spatial, Metabolomics 以及 Imaging 数据。",
        "lbl_language": "语言 / Language",

        # Query Modes
        "qmode_label": "🔍 选择查询模式",
        "qmode_gene": "🧬 全库搜索 (By Gene)",
        "qmode_scrna": "🔬 scRNA 细胞类型分析",
        "qmode_atac": "🧬 ATAC 样本开放度分析",
        "qmode_metabo": "⚗️ 代谢物关联分析",
        "qmode_spatial": "🗺️ 空间区域异质性分析",

        # Admin Tabs
        "tab_user_mgmt": "👥 用户管理",
        "tab_audit_logs": "📝 审计日志",
        "tab_data_crud": "🧬 组学数据维护",
        "tab_backup": "💾 备份与恢复",

        # User Mgmt
        "mgmt_create_user": "创建新用户",
        "mgmt_all_users": "用户列表",
        "lbl_username": "用户名",
        "lbl_password": "密码",
        "lbl_role": "权限角色",
        "btn_create_user": "创建用户",

        # CRUD General
        "crud_exp_anno": "📝 专家注释管理 (MySQL)",
        "crud_header_core": "🛠️ 核心组学数据修正 (SQLite)",
        "crud_select_db": "选择目标数据库",
        "crud_op_create": "➕ 新增 (Create)",
        "crud_op_update": "📝 修改 (Update)",
        "crud_op_delete": "🗑️ 删除 (Delete)",
        "btn_submit": "提交",
        "btn_add": "添加记录",
        "btn_update": "确认修改",
        "btn_delete": "确认删除",
        "msg_success": "操作成功！",

        # CRUD Fields (Specifics)
        "lbl_target_gene": "目标基因 (Gene)",
        "lbl_content": "注释内容",
        "lbl_subtype": "亚型 (Subtype)",
        "lbl_celltype": "细胞类型 (CellType)",
        "lbl_exp_val": "表达量 (Value)",
        "lbl_sample_id": "样本ID (Sample)",
        "lbl_region": "空间区域 (Region)",
        "lbl_metabo": "代谢物名称",
        "lbl_json": "JSON 内容",
        "lbl_old_val": "原值 (Old)",
        "lbl_new_val": "新值 (New)",

        # Backup
        "backup_title": "📦 全系统备份 (Full Backup)",
        "backup_desc": "将用户数据(MySQL)与组学数据(SQLite)打包下载。",
        "backup_sel_content": "1. 选择备份内容",
        "backup_inc_mysql": "用户数据 (MySQL)",
        "backup_inc_omics": "组学数据 (SQLite)",
        "btn_gen_zip": "📦 生成备份包",
        "restore_title": "♻️ 灾难恢复 (Disaster Restore)",
        "restore_sec": "2. 系统恢复",
        "restore_warn": "⚠️ 警告：上传的 ZIP 包将覆盖现有数据库文件，此操作不可逆！",
        "btn_start_restore": "🔥 开始恢复",

        # Query UI General
        "search_label": "输入基因 Symbol",
        "filter_label": "亚型过滤",
        "warn_no_data": "未找到相关数据。",
        "success_login": "验证通过！正在进入系统...",
        "err_login": "用户名或密码错误",
        "btn_pdf": "📄 生成分析报告 (Export PDF)",
        "info_expert_anno": "📋 专家注释 (Expert Annotations)",

        # Tabs
        "tab_scrna": "🔬 scRNA (单细胞)",
        "tab_atac": "🧬 ATAC (表观)",
        "tab_metabo": "⚗️ Metabo (代谢)",
        "tab_spatial": "🗺️ Spatial (空间)",
        "tab_imaging": "🖼️ Imaging (影像)",

        "data_browser": "📚 数据字典导览",
        "top_genes_list": "🔥 高表达基因",
        "top_metas_list": "🧪 高表达代谢物",
        "atac_sim_note": "⚠️ 注：当前 ATAC 数据库缺失临床亚型标注。下图展示基于模拟元数据的分组对比。",
        "atac_raw_title": "2. 原始样本分布 (未过滤)",
        "spatial_single_note": "ℹ️ 提示：当前 Spatial 模块展示标准参考样本 V1 (HER2_Positive)。",
        "imaging_note": "💡 说明：展示 AI 辅助识别的肿瘤感兴趣区域 (ROI)。",
        "input_gene_ph": "尝试: FOXA1, ESR1, PKM",
        "input_top_n": "筛选 Top N 结果",
        "header_top_genes": "🔥 高表达基因排行",
        "caption_plot_limit": "注：为保证图表清晰，图表仅展示前 50 项，完整数据请见下方表格。"
    },
    "EN": {
        "title": "BC-MOD Multi-Omics Database",
        "sidebar_header_user": "👤 User Center",
        "sidebar_header_nav": "🧭 Navigation",
        "sidebar_header_settings": "⚙️ Settings & Tools",
        "nav_query": "🔍 Data Query",
        "nav_admin": "🛠️ Admin Panel",
        "login_title": "🔐 System Login",
        "login_btn": "Login",
        "logout_btn": "Logout",
        "intro_title": "Introduction",
        "intro_text": "A database integrating scRNA-seq, ATAC-seq, Spatial, Metabolomics and Imaging data.",
        "lbl_language": "Language",

        "qmode_label": "🔍 Query Mode",
        "qmode_gene": "🧬 Global Search (By Gene)",
        "qmode_scrna": "🔬 scRNA Cell Type Analysis",
        "qmode_atac": "🧬 ATAC Sample Analysis",
        "qmode_metabo": "⚗️ Metabolite Analysis",
        "qmode_spatial": "🗺️ Spatial Region Analysis",

        "tab_user_mgmt": "👥 User Mgmt",
        "tab_audit_logs": "📝 Audit Logs",
        "tab_data_crud": "🧬 Data Maintenance",
        "tab_backup": "💾 Backup & Restore",

        "mgmt_create_user": "Create New User",
        "mgmt_all_users": "All Users",
        "lbl_username": "Username",
        "lbl_password": "Password",
        "lbl_role": "Role",
        "btn_create_user": "Create User",

        "crud_exp_anno": "📝 Expert Annotations (MySQL)",
        "crud_header_core": "🛠️ Core Omics Maintenance (SQLite)",
        "crud_select_db": "Select Database",
        "crud_op_create": "➕ Create",
        "crud_op_update": "📝 Update",
        "crud_op_delete": "🗑️ Delete",
        "btn_submit": "Submit",
        "btn_add": "Add Record",
        "btn_update": "Update",
        "btn_delete": "Delete",
        "msg_success": "Operation Successful!",

        "lbl_target_gene": "Target Gene",
        "lbl_content": "Content",
        "lbl_subtype": "Subtype",
        "lbl_celltype": "CellType",
        "lbl_exp_val": "Expression Value",
        "lbl_sample_id": "Sample ID",
        "lbl_region": "Region",
        "lbl_metabo": "Metabolite Name",
        "lbl_json": "JSON Content",
        "lbl_old_val": "Old Value",
        "lbl_new_val": "New Value",

        "backup_title": "📦 Full System Backup",
        "backup_desc": "Download User Data (MySQL) and Omics Data (SQLite) as ZIP.",
        "backup_sel_content": "1. Select Content",
        "backup_inc_mysql": "User Data (MySQL)",
        "backup_inc_omics": "Omics Data (SQLite)",
        "btn_gen_zip": "📦 Generate ZIP",
        "restore_title": "♻️ Disaster Restore",
        "restore_sec": "2. System Restore",
        "restore_warn": "⚠️ Warning: This will overwrite existing databases!",
        "btn_start_restore": "🔥 Start Restore",

        "search_label": "Enter Gene Symbol",
        "filter_label": "Subtype Filter",
        "warn_no_data": "No data found.",
        "success_login": "Login successful! Redirecting...",
        "err_login": "Invalid username or password",
        "btn_pdf": "📄 Export PDF Report",
        "info_expert_anno": "📋 Expert Annotations",

        "tab_scrna": "🔬 scRNA",
        "tab_atac": "🧬 ATAC",
        "tab_metabo": "⚗️ Metabo",
        "tab_spatial": "🗺️ Spatial",
        "tab_imaging": "🖼️ Imaging",

        "data_browser": "📚 Data Dictionary",
        "top_genes_list": "🔥 Top Genes (scRNA)",
        "top_metas_list": "🧪 Top Metabolites",
        "atac_sim_note": "⚠️ Note: ATAC subtypes are simulated for demonstration.",
        "atac_raw_title": "2. Raw Sample Distribution",
        "spatial_single_note": "ℹ️ Note: Showing Reference Sample V1 (HER2_Positive).",
        "imaging_note": "💡 Note: AI-identified Tumor ROI demo.",
        "input_gene_ph": "Try: FOXA1, ESR1, PKM",
        "input_top_n": "Show Top N",
        "header_top_genes": "🔥 Top Expressed Genes",
        "caption_plot_limit": "Note: Plot limited to top 50 items for clarity. See table for full list."
    }
}

# ==========================================
# 4. 配置与工具函数
# ==========================================
st.set_page_config(page_title="BC-MOD Database", page_icon="🧬", layout="wide")

DB_PATHS = {
    "scRNA": "./dbs/scrna_3nf.db",
    "ATAC": "./dbs/atac.db",
    "Metabo": "./dbs/metabolomics.db",
    "Spatial": "./dbs/spatial.db",
    "Imaging": "./dbs/imaging.db"
}

if 'logged_in' not in st.session_state:
    st.session_state['logged_in'] = False
    st.session_state['user_role'] = None
    st.session_state['username'] = None
if 'lang' not in st.session_state:
    st.session_state['lang'] = "CN"


def t(key):
    """Retrieve translation safely."""
    return TRANS[st.session_state['lang']].get(key, key)


def run_sqlite_query(db_key, sql):
    db_path = DB_PATHS.get(db_key)
    if not os.path.exists(db_path): return None
    try:
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(sql, conn)
        conn.close()
        return df
    except:
        return None


def get_atac_meta(sample_id):
    types = ["TNBC", "HER2_Positive", "Luminal_A", "Luminal_B", "Normal"]
    return types[hash(sample_id) % len(types)]


def get_distinct_values(db_key, table, column):
    sql = f"SELECT DISTINCT {column} FROM {table} ORDER BY {column}"
    df = run_sqlite_query(db_key, sql)
    return df[column].tolist() if df is not None else []


def get_real_top_elements(mode):
    try:
        if mode == "gene":
            sql = "SELECT Gene FROM Table_Expression GROUP BY Gene ORDER BY SUM(Avg_Expression) DESC LIMIT 10"
            df = run_sqlite_query("scRNA", sql)
            return df['Gene'].tolist() if df is not None else []
        elif mode == "metabo":
            sql = "SELECT Metabolite FROM Metabolite_Expression GROUP BY Metabolite ORDER BY SUM(Expression_Level) DESC LIMIT 10"
            df = run_sqlite_query("Metabo", sql)
            return df['Metabolite'].tolist() if df is not None else []
    except:
        return []
    return []


def ensure_annotation_table():
    conn = auth.get_connection()
    if conn:
        cursor = conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS gene_annotations (
                id INT AUTO_INCREMENT PRIMARY KEY,
                gene VARCHAR(50),
                note TEXT,
                author VARCHAR(50),
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        conn.commit()
        conn.close()


def safe_cursor_fetch(conn, sql, params=None):
    try:
        cursor = conn.cursor()
        if params:
            cursor.execute(sql, params)
        else:
            cursor.execute(sql)
        if cursor.description:
            cols = [col[0] for col in cursor.description]
            data = cursor.fetchall()
            return pd.DataFrame(data, columns=cols)
        return pd.DataFrame()
    except:
        return pd.DataFrame()


# --- Enhanced PDF Generation ---
def create_pdf_report(gene_name, username, query_data):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(0, 10, txt="BC-MOD Multi-Omics Analysis Report", ln=1, align='C')

    pdf.set_font("Arial", size=10)
    pdf.cell(0, 10, txt=f"Generated by: {username} | Date: {time.strftime('%Y-%m-%d %H:%M')}", ln=1, align='C')
    pdf.line(10, 30, 200, 30)
    pdf.ln(10)

    pdf.set_font("Arial", 'B', 12)
    pdf.cell(0, 10, txt=f"Target Gene: {gene_name}", ln=1)

    pdf.set_font("Arial", size=11)
    for omics, summary in query_data.items():
        pdf.set_font("Arial", 'B', 11)
        pdf.set_fill_color(240, 240, 240)
        pdf.cell(0, 8, txt=f"[{omics}] Analysis Results:", ln=1, fill=True)
        pdf.set_font("Arial", size=10)
        if summary:
            pdf.multi_cell(0, 6, txt=summary)
        else:
            pdf.multi_cell(0, 6, txt="No data available.")
        pdf.ln(4)

    pdf.ln(10)
    pdf.set_font("Arial", 'I', 8)
    pdf.cell(0, 10, txt="BC-MOD System | Database Course Design 2025", align='C')
    return pdf.output(dest='S').encode('latin-1')


# ==========================================
# 5. UI Modules
# ==========================================

def login_ui():
    st.markdown(f"## 🧬 {t('title')}")
    st.markdown("---")
    lang = st.radio("Language / 语言", ["CN", "EN"], horizontal=True)
    st.session_state['lang'] = lang

    c1, c2 = st.columns([1, 1])
    with c1:
        st.info(f"**{t('intro_text')}**\n\n- Admin: `admin` / `admin123456`\n- Guest: `guest` / `guest123456`")
    with c2:
        st.subheader(t('login_title'))
        user = st.text_input(t('lbl_username'))
        pwd = st.text_input(t('lbl_password'), type="password")
        if st.button(t('login_btn'), type="primary", use_container_width=True):
            u = auth.check_login(user, pwd)
            if u:
                st.session_state['logged_in'] = True
                st.session_state['user_role'] = u['role']
                st.session_state['username'] = u['username']
                auth.log_action(user, "Login")
                ensure_annotation_table()
                st.success(t('success_login'))
                st.rerun()
            else:
                st.error(t('err_login'))


def admin_ui():
    st.header(t('nav_admin'))
    st.warning(f"Admin: {st.session_state['username']}")

    tab1, tab2, tab3, tab4 = st.tabs([t('tab_user_mgmt'), t('tab_audit_logs'), t('tab_data_crud'), t('tab_backup')])

    # --- Tab 1: User Management ---
    with tab1:
        c1, c2 = st.columns([1, 2])
        with c1:
            st.subheader(t('mgmt_create_user'))
            with st.form("new_u"):
                nu = st.text_input(t('lbl_username'))
                np = st.text_input(t('lbl_password'), type="password")
                nr = st.selectbox(t('lbl_role'), ["guest", "admin"])
                if st.form_submit_button(t('btn_create_user')):
                    if auth.create_user(nu, np, nr):
                        st.success(f"User {nu} created!")
                        auth.log_action(st.session_state['username'], f"Create user {nu}")
                        st.rerun()
                    else:
                        st.error("Failed.")

        with c2:
            st.subheader(t('mgmt_all_users'))
            conn = auth.get_connection()
            if conn:
                df_users = safe_cursor_fetch(conn, "SELECT id, username, role, created_at FROM users")
                conn.close()
                st.dataframe(df_users, use_container_width=True, hide_index=True)

    # --- Tab 2: Logs ---
    with tab2:
        if st.button("Refresh Logs"): st.rerun()
        st.dataframe(auth.get_system_logs(50), use_container_width=True)

    # --- Tab 3: Unified CRUD ---
    with tab3:
        st.subheader(t('tab_data_crud'))

        # 1. Expert Annotations
        with st.expander(t('crud_exp_anno'), expanded=True):
            c1, c2 = st.columns([1, 2])
            target_gene = c1.text_input(t('lbl_target_gene'), "FOXA1").strip().upper()
            note_content = c2.text_input(t('lbl_content'))
            if st.button(t('btn_submit')):
                conn = auth.get_connection()
                if conn:
                    cursor = conn.cursor()
                    cursor.execute("INSERT INTO gene_annotations (gene, note, author) VALUES (%s, %s, %s)",
                                   (target_gene, note_content, st.session_state['username']))
                    conn.commit()
                    conn.close()
                    st.success(t('msg_success'))
                    auth.log_action(st.session_state['username'], f"Create Annotation: {target_gene}")

        st.divider()

        # 2. Omics Data Maintenance
        st.markdown(f"#### {t('crud_header_core')}")
        db_options = ["scRNA", "ATAC", "Metabo", "Spatial", "Imaging"]
        selected_db = st.selectbox(t('crud_select_db'), db_options)

        with st.container(border=True):
            st.markdown(f"**{selected_db} Maintenance**")
            crud_tab1, crud_tab2, crud_tab3 = st.tabs([t('crud_op_create'), t('crud_op_update'), t('crud_op_delete')])

            # --- scRNA CRUD ---
            if selected_db == "scRNA":
                conn = sqlite3.connect(DB_PATHS['scRNA'])
                with crud_tab1:  # Create
                    c1, c2, c3, c4 = st.columns(4)
                    new_gene = c1.text_input(t('lbl_target_gene'), key="c_sc_g").upper()
                    new_sub = c2.selectbox(t('lbl_subtype'),
                                           ["TNBC", "Luminal_A", "HER2_Positive", "Normal", "Luminal_B"], key="c_sc_s")
                    new_cell = c3.text_input(t('lbl_celltype'), "T_cells", key="c_sc_c")
                    new_val = c4.number_input(t('lbl_exp_val'), 0.0, 1000.0, 1.0, key="c_sc_v")
                    if st.button(t('btn_add'), key="scrna_add"):
                        try:
                            cur = conn.cursor()
                            cur.execute("SELECT gene_id FROM Genes WHERE gene_name=?", (new_gene,))
                            res = cur.fetchone()
                            gid = res[0] if res else cur.execute("INSERT INTO Genes (gene_name) VALUES (?)",
                                                                 (new_gene,)).lastrowid

                            cur.execute("SELECT group_id FROM CellGroups WHERE subtype=? AND celltype=?",
                                        (new_sub, new_cell))
                            res = cur.fetchone()
                            grid = res[0] if res else cur.execute(
                                "INSERT INTO CellGroups (subtype, celltype) VALUES (?,?)",
                                (new_sub, new_cell)).lastrowid

                            cur.execute("INSERT INTO Expression (gene_id, group_id, value) VALUES (?,?,?)",
                                        (gid, grid, new_val))
                            conn.commit()
                            st.success(t('msg_success'))
                            auth.log_action(st.session_state['username'], f"Create scRNA: {new_gene}")
                        except Exception as e:
                            st.error(f"Error: {e}")

                with crud_tab2:  # Update
                    c1, c2, c3, c4 = st.columns(4)
                    u_gene = c1.text_input(t('lbl_target_gene'), "FOXA1", key="u_sc_g").upper()
                    u_sub = c2.selectbox(t('lbl_subtype'), ["TNBC", "Luminal_A"], key="u_sc_s")
                    u_cell = c3.text_input(t('lbl_celltype'), "T_cells", key="u_sc_c")
                    u_val = c4.number_input(t('lbl_new_val'), 0.0, 1000.0, 5.0, key="u_sc_v")
                    if st.button(t('btn_update'), key="scrna_upd"):
                        try:
                            cur = conn.cursor()
                            cur.execute("SELECT gene_id FROM Genes WHERE gene_name=?", (u_gene,))
                            res_g = cur.fetchone()
                            cur.execute("SELECT group_id FROM CellGroups WHERE subtype=? AND celltype=?",
                                        (u_sub, u_cell))
                            res_grp = cur.fetchone()
                            if res_g and res_grp:
                                cur.execute("UPDATE Expression SET value=? WHERE gene_id=? AND group_id=?",
                                            (u_val, res_g[0], res_grp[0]))
                                conn.commit()
                                st.success(t('msg_success'))
                            else:
                                st.warning("Record not found.")
                        except Exception as e:
                            st.error(f"{e}")

                with crud_tab3:  # Delete
                    del_gene = st.text_input("Delete Gene (Symbol)", key="scrna_del").upper()
                    if st.button(t('btn_delete'), key="scrna_del_btn"):
                        try:
                            cur = conn.cursor()
                            cur.execute("SELECT gene_id FROM Genes WHERE gene_name=?", (del_gene,))
                            res = cur.fetchone()
                            if res:
                                gid = res[0]
                                cur.execute("DELETE FROM Expression WHERE gene_id=?", (gid,))
                                cur.execute("DELETE FROM Genes WHERE gene_id=?", (gid,))
                                conn.commit()
                                st.success(t('msg_success'))
                            else:
                                st.warning("Not found.")
                        except Exception as e:
                            st.error(f"{e}")
                conn.close()

            # --- ATAC CRUD ---
            elif selected_db == "ATAC":
                conn = sqlite3.connect(DB_PATHS['ATAC'])
                with crud_tab1:
                    c1, c2 = st.columns(2)
                    new_samp = c1.text_input(t('lbl_sample_id'), key="at_c_s")
                    ref_gene = c2.text_input("Ref Gene", "FOXA1", key="at_c_g")
                    val = st.number_input("Value", 0.0, key="at_c_v")
                    if st.button(t('btn_add'), key="atac_add"):
                        try:
                            cur = conn.cursor()
                            cur.execute(f"INSERT INTO sample_gene_matrix (sample, {ref_gene}) VALUES (?, ?)",
                                        (new_samp, val))
                            conn.commit()
                            st.success(t('msg_success'))
                        except Exception as e:
                            st.error(f"{e}")

                with crud_tab2:
                    c1, c2, c3 = st.columns(3)
                    samples = pd.read_sql("SELECT sample FROM sample_gene_matrix", conn)['sample'].tolist()
                    tgt_sample = c1.selectbox(t('lbl_sample_id'), samples, key="at_u_s")
                    tgt_gene = c2.text_input("Gene", "FOXA1", key="at_u_g").upper()
                    new_val = c3.number_input(t('lbl_new_val'), 0.0, key="at_u_v")
                    if st.button(t('btn_update'), key="atac_upd"):
                        try:
                            cursor = conn.cursor()
                            cursor.execute(f"UPDATE sample_gene_matrix SET {tgt_gene} = ? WHERE sample = ?",
                                           (new_val, tgt_sample))
                            conn.commit()
                            st.success(t('msg_success'))
                        except Exception as e:
                            st.error(f"{e}")

                with crud_tab3:
                    del_sample = st.selectbox(t('lbl_sample_id'), samples, key="atac_del")
                    if st.button(t('btn_delete'), key="atac_del_btn"):
                        conn.execute("DELETE FROM sample_gene_matrix WHERE sample = ?", (del_sample,))
                        conn.commit()
                        st.success(t('msg_success'))
                conn.close()

            # --- Metabo CRUD ---
            elif selected_db == "Metabo":
                conn = sqlite3.connect(DB_PATHS['Metabo'])
                with crud_tab1:
                    c1, c2, c3 = st.columns(3)
                    m_name = c1.text_input(t('lbl_metabo'), key="mt_c_n")
                    m_sub = c2.selectbox(t('lbl_subtype'), ["TNBC", "Normal"], key="mt_c_s")
                    m_val = c3.number_input(t('lbl_exp_val'), key="mt_c_v")
                    if st.button(t('btn_add'), key="met_add"):
                        conn.execute("INSERT INTO Metabolite_Expression VALUES (?,?,?)", (m_name, m_sub, m_val))
                        conn.commit()
                        st.success(t('msg_success'))

                with crud_tab2:
                    c1, c2, c3 = st.columns(3)
                    metas = pd.read_sql("SELECT DISTINCT Metabolite FROM Metabolite_Expression", conn)[
                        'Metabolite'].tolist()
                    u_meta = c1.selectbox(t('lbl_metabo'), metas, key="mt_u_n")
                    u_sub = c2.selectbox(t('lbl_subtype'), ["TNBC", "Normal"], key="mt_u_s")
                    u_val = c3.number_input(t('lbl_new_val'), key="mt_u_v")
                    if st.button(t('btn_update'), key="met_upd"):
                        conn.execute(
                            "UPDATE Metabolite_Expression SET Expression_Level=? WHERE Metabolite=? AND Subtype=?",
                            (u_val, u_meta, u_sub))
                        conn.commit()
                        st.success(t('msg_success'))

                with crud_tab3:
                    del_meta = st.selectbox(t('lbl_metabo'), metas, key="met_del")
                    if st.button(t('btn_delete'), key="met_del_btn"):
                        conn.execute("DELETE FROM Metabolite_Expression WHERE Metabolite = ?", (del_meta,))
                        conn.commit()
                        st.success(t('msg_success'))
                conn.close()

            # --- Spatial CRUD ---
            elif selected_db == "Spatial":
                conn = sqlite3.connect(DB_PATHS['Spatial'])
                with crud_tab1:
                    c1, c2, c3, c4 = st.columns(4)
                    g = c1.text_input(t('lbl_target_gene'), key="sp_c_g").upper()
                    s = c2.text_input(t('lbl_sample_id'), "V1_...", key="sp_c_s")
                    r = c3.text_input(t('lbl_region'), "Region 1", key="sp_c_r")
                    v = c4.number_input(t('lbl_exp_val'), key="sp_c_v")
                    if st.button(t('btn_add'), key="sp_add"):
                        conn.execute("INSERT INTO Table_SpatialExpression VALUES (?,?,?,?)", (g, s, r, v))
                        conn.commit()
                        st.success(t('msg_success'))

                with crud_tab2:
                    c1, c2, c3 = st.columns(3)
                    u_g = c1.text_input(t('lbl_target_gene'), "FOXA1", key="sp_u_g").upper()
                    u_r = c2.text_input(t('lbl_region'), "Region 1", key="sp_u_r")
                    u_v = c3.number_input(t('lbl_new_val'), key="sp_u_v")
                    if st.button(t('btn_update'), key="sp_upd"):
                        conn.execute("UPDATE Table_SpatialExpression SET Avg_Expression=? WHERE Gene=? AND Region=?",
                                     (u_v, u_g, u_r))
                        conn.commit()
                        st.success(t('msg_success'))

                with crud_tab3:
                    d_g = st.text_input(t('lbl_target_gene'), key="sp_del").upper()
                    if st.button(t('btn_delete'), key="sp_del_btn"):
                        conn.execute("DELETE FROM Table_SpatialExpression WHERE Gene=?", (d_g,))
                        conn.commit()
                        st.success(t('msg_success'))
                conn.close()

            # --- Imaging CRUD ---
            elif selected_db == "Imaging":
                conn = sqlite3.connect(DB_PATHS['Imaging'])
                with crud_tab1:
                    st.caption("Insert raw JSON")
                    json_str = st.text_area(t('lbl_json'), '{"positive": []}')
                    if st.button(t('btn_add'), key="img_add"):
                        try:
                            json.loads(json_str)
                            conn.execute("INSERT INTO annotations (annotation) VALUES (?)", (json_str,))
                            conn.commit()
                            st.success(t('msg_success'))
                        except:
                            st.error("Invalid JSON")

                with crud_tab2:
                    cursor = conn.cursor()
                    cursor.execute("SELECT id, annotation FROM annotations ORDER BY id DESC LIMIT 1")
                    res = cursor.fetchone()
                    if res:
                        old_id, old_json = res
                        new_json = st.text_area(f"Edit ID {old_id}", old_json, height=150)
                        if st.button(t('btn_update'), key="img_upd"):
                            try:
                                json.loads(new_json)
                                conn.execute("UPDATE annotations SET annotation=? WHERE id=?", (new_json, old_id))
                                conn.commit()
                                st.success(t('msg_success'))
                            except:
                                st.error("Invalid JSON")
                    else:
                        st.info("No annotations.")

                with crud_tab3:
                    id_to_del = st.number_input("Annotation ID", 1, step=1, key="img_del_id")
                    if st.button(t('btn_delete'), key="img_del"):
                        conn.execute("DELETE FROM annotations WHERE id=?", (id_to_del,))
                        conn.commit()
                        st.success(t('msg_success'))
                conn.close()

    # --- Tab 4: Backup ---
    with tab4:
        st.subheader(t('backup_title'))
        st.caption(t('backup_desc'))

        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"##### {t('backup_sel_content')}")
            inc_mysql = st.checkbox(t('backup_inc_mysql'), value=True)
            inc_omics = st.multiselect(t('backup_inc_omics'), options=list(DB_PATHS.keys()), default=["scRNA"])

            if st.button(t('btn_gen_zip')):
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                    if inc_mysql:
                        conn = auth.get_connection()
                        if conn:
                            for table in ['users', 'system_logs', 'gene_annotations']:
                                df = safe_cursor_fetch(conn, f"SELECT * FROM {table}")
                                if not df.empty:
                                    zip_file.writestr(f"mysql_backup/{table}.csv", df.to_csv(index=False))
                            conn.close()
                    for db_key in inc_omics:
                        path = DB_PATHS.get(db_key)
                        if os.path.exists(path):
                            zip_file.write(path, arcname=f"sqlite_backup/{os.path.basename(path)}")

                st.download_button("📥 Download", data=zip_buffer.getvalue(), file_name="bc_mod_backup.zip",
                                   mime="application/zip")

        with c2:
            st.markdown(f"##### {t('restore_sec')}")
            st.warning(t('restore_warn'))
            uploaded_zip = st.file_uploader("ZIP File", type="zip")
            if uploaded_zip and st.button(t('btn_start_restore')):
                try:
                    with zipfile.ZipFile(uploaded_zip, 'r') as z:
                        for file in z.namelist():
                            if file.startswith("sqlite_backup/") and file.endswith(".db"):
                                filename = os.path.basename(file)
                                target_path = None
                                for k, v in DB_PATHS.items():
                                    if os.path.basename(v) == filename:
                                        target_path = v
                                        break
                                if target_path:
                                    with open(target_path, "wb") as f:
                                        f.write(z.read(file))
                                    st.info(f"Restored: {filename}")
                    st.success(t('msg_success'))
                except Exception as e:
                    st.error(f"Error: {e}")


# ------------------------------------------
# Main Query Logic
# ------------------------------------------
def query_ui():
    with st.sidebar:
        st.markdown(f"### {t('sidebar_header_user')}")
        st.success(f"**{st.session_state['username']}** ({st.session_state['user_role']})")
        if st.button(t('logout_btn'), use_container_width=True):
            auth.log_action(st.session_state['username'], "Logout")
            st.session_state['logged_in'] = False
            st.rerun()

        st.markdown("---")
        st.markdown(f"### {t('sidebar_header_nav')}")
        nav_mode = st.radio("Go to", [t('nav_query'), t('nav_admin')], label_visibility="collapsed")

        if nav_mode == t('nav_admin'):
            if st.session_state['user_role'] == 'admin':
                return admin_ui()
            else:
                st.error("Access Denied")
                return

        st.markdown("---")
        with st.expander(t('sidebar_header_settings'), expanded=False):
            current_idx = 0 if st.session_state['lang'] == "CN" else 1
            new_lang = st.selectbox(t('lbl_language'), ["CN", "EN"], index=current_idx)
            if new_lang != st.session_state['lang']:
                st.session_state['lang'] = new_lang
                st.rerun()

        with st.expander(t('data_browser'), expanded=False):
            st.markdown(f"**{t('top_genes_list')}**")
            genes = get_real_top_elements("gene")
            st.code("\n".join(genes) if genes else "Loading...")
            st.markdown(f"**{t('top_metas_list')}**")
            metas = get_real_top_elements("metabo")
            st.code("\n".join(metas) if metas else "Loading...")

    st.title(f"🧬 {t('title')}")

    q_modes = {
        t('qmode_gene'): "gene",
        t('qmode_scrna'): "scrna_advanced",
        t('qmode_atac'): "atac_advanced",
        t('qmode_metabo'): "metabo_advanced",
        t('qmode_spatial'): "spatial_advanced"
    }

    col_mode, col_blank = st.columns([1, 2])
    with col_mode:
        selected_mode_label = st.selectbox(t('qmode_label'), list(q_modes.keys()))

    mode_key = q_modes[selected_mode_label]
    pdf_summary_data = {}

    if mode_key == "gene":
        c1, c2 = st.columns([3, 1])
        gene_input = c1.text_input(t('search_label'), value="FOXA1", placeholder=t('input_gene_ph')).strip().upper()
        subtype = c2.selectbox(t('filter_label'), ["All", "TNBC", "HER2_Positive", "Luminal_A", "Luminal_B", "Normal"])

        if not gene_input: return

        tabs = st.tabs([t('tab_scrna'), t('tab_atac'), t('tab_metabo'), t('tab_spatial'), t('tab_imaging')])

        with tabs[0]:
            conn_mysql = auth.get_connection()
            if conn_mysql:
                df_notes = safe_cursor_fetch(conn_mysql, "SELECT * FROM gene_annotations WHERE gene = %s",
                                             (gene_input,))
                conn_mysql.close()
                if not df_notes.empty:
                    st.info(f"{t('info_expert_anno')} ({len(df_notes)})")
                    for idx, row in df_notes.iterrows():
                        st.markdown(f"- {row['note']} *(By: {row['author']})*")

            sql = f"SELECT Subtype, CellType, Avg_Expression FROM Table_Expression WHERE Gene = '{gene_input}'"
            if subtype != "All": sql += f" AND Subtype = '{subtype}'"
            df = run_sqlite_query("scRNA", sql)
            if df is not None and not df.empty:
                st.bar_chart(df, x="CellType", y="Avg_Expression", color="Subtype")
                pdf_summary_data['scRNA'] = f"Found {len(df)} records.\nMax: {df['Avg_Expression'].max():.2f}"
            else:
                st.warning(t('warn_no_data'))

        with tabs[1]:
            sql = f"SELECT sample, {gene_input} FROM sample_gene_matrix"
            df = run_sqlite_query("ATAC", sql)
            if df is not None and not df.empty:
                df['Subtype'] = df['sample'].apply(get_atac_meta)
                st.info(t('atac_sim_note'))
                df_filter = df[df['Subtype'] == subtype] if subtype != "All" else df
                avg = df_filter.groupby("Subtype")[gene_input].mean().reset_index()
                st.bar_chart(avg, x="Subtype", y=gene_input, color="Subtype")
                st.markdown("---")
                st.bar_chart(df, x="sample", y=gene_input)
                pdf_summary_data['ATAC'] = f"Accession samples: {len(df)}\nAvg Openness: {df[gene_input].mean():.2f}"
            else:
                st.info(f"Gene {gene_input} not found in ATAC DB.")

        with tabs[2]:
            df_map = run_sqlite_query("Metabo", f"SELECT * FROM Gene_Metabolite_Map WHERE Gene = '{gene_input}'")
            if df_map is not None and not df_map.empty:
                st.dataframe(df_map)
                metas = [m.replace("'", "''") for m in df_map['Metabolite'].tolist()]
                if metas:
                    m_str = "', '".join(metas)
                    sql = f"SELECT * FROM Metabolite_Expression WHERE Metabolite IN ('{m_str}')"
                    if subtype != "All": sql += f" AND Subtype = '{subtype}'"
                    df_exp = run_sqlite_query("Metabo", sql)
                    if df_exp is not None and not df_exp.empty:
                        st.line_chart(df_exp, x="Subtype", y="Expression_Level", color="Metabolite")
            else:
                st.warning(t('warn_no_data'))

        with tabs[3]:
            st.info(t('spatial_single_note'))
            sql = f"SELECT SampleID, Region, Avg_Expression FROM Table_SpatialExpression WHERE Gene = '{gene_input}'"
            df = run_sqlite_query("Spatial", sql)
            if df is not None and not df.empty:
                st.bar_chart(df, x="Region", y="Avg_Expression", color="SampleID")
                pdf_summary_data['Spatial'] = "Sample V1 Detected."
            else:
                st.warning(t('warn_no_data'))

        with tabs[4]:
            st.info(t('imaging_note'))
            df_anno = run_sqlite_query("Imaging", "SELECT annotation FROM annotations LIMIT 1")
            if df_anno is not None:
                try:
                    data = json.loads(df_anno.iloc[0]['annotation'])
                    fig, ax = plt.subplots(figsize=(6, 6))
                    ax.set_title("Pathology ROI Annotation")
                    ax.set_xlim(5000, 22000);
                    ax.set_ylim(13000, 8000)
                    for region in data.get('positive', []):
                        poly = patches.Polygon(region['vertices'], closed=True, facecolor='#FF4B4B', alpha=0.4,
                                               edgecolor='red')
                        ax.add_patch(poly)
                    st.pyplot(fig)
                except:
                    st.error("JSON Error")
            else:
                st.info("No Imaging Data")

        st.markdown("---")
        if st.button(t('btn_pdf')):
            pdf_bytes = create_pdf_report(gene_input, st.session_state['username'], pdf_summary_data)
            b64 = base64.b64encode(pdf_bytes).decode()
            href = f'<a href="data:application/octet-stream;base64,{b64}" download="BC_MOD_{gene_input}_Report.pdf">📥 Download PDF</a>'
            st.markdown(href, unsafe_allow_html=True)
            auth.log_action(st.session_state['username'], f"Export PDF: {gene_input}")

    # --- Other Modes (Shortened) ---
    elif mode_key == "scrna_advanced":
        st.subheader(t('qmode_scrna'))
        cell_types = get_distinct_values("scRNA", "Table_Expression", "CellType")
        c1, c2 = st.columns(2)
        ct = c1.selectbox(t('input_celltype'), cell_types)
        top_n = c2.slider(t('input_top_n'), 3, 100, 5)
        if ct:
            sql = f"SELECT Gene, AVG(Avg_Expression) as MeanExpr FROM Table_Expression WHERE CellType = '{ct}' GROUP BY Gene ORDER BY MeanExpr DESC LIMIT {top_n}"
            df = run_sqlite_query("scRNA", sql)
            if df is not None:
                st.bar_chart(df.head(50), x="Gene", y="MeanExpr")
                st.dataframe(df, use_container_width=True)

    elif mode_key == "atac_advanced":
        st.subheader(t('qmode_atac'))
        samples = get_distinct_values("ATAC", "sample_gene_matrix", "sample")
        sid = st.selectbox(t('input_sample'), samples)
        if sid:
            df = run_sqlite_query("ATAC", f"SELECT * FROM sample_gene_matrix WHERE sample = '{sid}'")
            if df is not None:
                df_t = df.drop(columns=['sample']).T
                st.bar_chart(df_t.head(50))

    elif mode_key == "metabo_advanced":
        st.subheader(t('qmode_metabo'))
        metas = get_distinct_values("Metabo", "Gene_Metabolite_Map", "Metabolite")
        m = st.selectbox(t('input_metabo'), metas)
        if m:
            m_safe = m.replace("'", "''")
            df = run_sqlite_query("Metabo",
                                  f"SELECT Subtype, Expression_Level FROM Metabolite_Expression WHERE Metabolite = '{m_safe}'")
            if df is not None: st.bar_chart(df, x="Subtype", y="Expression_Level")

    elif mode_key == "spatial_advanced":
        st.subheader(t('qmode_spatial'))
        regions = get_distinct_values("Spatial", "Table_SpatialExpression", "Region")
        r = st.selectbox(t('input_region'), regions)
        if r:
            df = run_sqlite_query("Spatial",
                                  f"SELECT Gene, Avg_Expression FROM Table_SpatialExpression WHERE Region = '{r}' ORDER BY Avg_Expression DESC LIMIT 10")
            if df is not None: st.bar_chart(df, x="Gene", y="Avg_Expression")


if __name__ == "__main__":
    if st.session_state['logged_in']:
        query_ui()
    else:
        login_ui()