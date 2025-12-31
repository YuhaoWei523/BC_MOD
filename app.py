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
# 3. 国际化字典 (i18n) - 修复与增强版
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
        "qmode_metabo": "⚗️ 代谢物高级分析",
        "qmode_spatial": "🗺️ 空间区域异质性分析",

        # Admin
        "tab_user_mgmt": "👥 用户管理",
        "tab_audit_logs": "📝 审计日志",
        "tab_data_crud": "🧬 组学数据维护",
        "tab_backup": "💾 备份与恢复",

        # CRUD
        "crud_exp_anno": "📝 专家注释管理 (MySQL)",
        "crud_header_core": "🛠️ 核心组学数据修正 (SQLite)",
        "crud_select_db": "选择目标数据库",
        "crud_mode_anno_add": "➕ 新增注释",
        "crud_mode_anno_manage": "🖊️ 管理已有注释 (修改/删除)",
        "crud_op_create": "➕ 新增数据 (Create)",
        "crud_op_update": "📝 修改/扩展 (Update/Extend)",
        "crud_op_delete": "🗑️ 删除数据 (Delete)",
        "crud_op_metabo_map": "🔗 基因-代谢物映射 (Mapping)",
        "crud_op_metabo_exp": "📊 代谢表达量 (Expression)",

        "btn_submit": "提交",
        "btn_save": "保存修改",
        "btn_delete": "删除记录",
        "btn_add": "执行添加",
        "btn_update": "执行更新",
        "btn_extend": "扩展列并更新",

        "msg_success": "操作成功！",
        "msg_fail": "操作失败，请检查输入或日志。",
        "msg_confirm_del": "确定要删除这条记录吗？",

        # Labels
        "lbl_target_gene": "目标基因 (Gene)",
        "lbl_content": "注释内容",
        "lbl_subtype": "亚型 (Subtype)",
        "lbl_celltype": "细胞类型 (CellType)",
        "lbl_exp_val": "数值 (Value)",
        "lbl_sample_id": "样本ID (Sample)",
        "lbl_ref_gene": "参考基因 (Ref Gene)",
        "lbl_col_gene": "基因列名 (Gene Column)",
        "lbl_region": "空间区域 (Region)",
        "lbl_metabo": "代谢物 (Metabolite)",
        "lbl_kegg": "KEGG 通路",
        "lbl_note": "备注 (Note)",
        "lbl_json": "JSON 内容 (ROI Data)",
        "lbl_anno_id": "注释ID",
        "lbl_new_sample": "新样本ID (New Sample)",
        "lbl_auto_col": "⚠️ 若基因列不存在，将自动修改表结构 (ALTER TABLE) 增加该列。",
        "lbl_metabo_mode": "维护模式",
        "lbl_search_kegg": "筛选 KEGG 通路",
        "lbl_json_example": "示例数据已加载，请基于此修改：",

        # Backup & Other
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

        "search_label": "输入基因 Symbol",
        "filter_label": "亚型过滤",
        "warn_no_data": "未找到相关数据。",
        "success_login": "验证通过！正在进入系统...",
        "err_login": "用户名或密码错误",
        "btn_pdf": "📄 生成分析报告 (Export PDF)",
        "info_expert_anno": "📋 专家注释 (Expert Annotations)",

        # Tabs & Notes (Fixing missed keys)
        "tab_scrna": "🔬 scRNA",
        "tab_atac": "🧬 ATAC",
        "tab_metabo": "⚗️ Metabo",
        "tab_spatial": "🗺️ Spatial",
        "tab_imaging": "🖼️ Imaging",
        "imaging_note": "💡 说明：展示 AI 辅助识别的肿瘤感兴趣区域 (ROI)。JSON 数据定义了多边形顶点坐标。",
        "spatial_single_note": "ℹ️ 提示：当前 Spatial 模块展示标准参考样本 V1 (HER2_Positive)。",
        "atac_sim_note": "⚠️ 注：ATAC 数据展示包含真实录入的亚型与模拟亚型。",
        "atac_raw_title": "2. 原始样本分布 (未过滤)",

        "data_browser": "📚 数据字典导览",
        "top_genes_list": "🔥 高表达基因 (scRNA)",
        "top_metas_list": "🧪 高表达代谢物",
        "input_gene_ph": "尝试: FOXA1, ESR1, PKM",
        "input_top_n": "显示 Top N",
        "input_sample": "选择样本",
        "input_metabo": "选择代谢物",
        "input_region": "选择区域",
        "input_celltype": "选择细胞类型",
        "header_top_genes": "🔥 高表达基因排行",
        "mgmt_create_user": "创建用户",
        "mgmt_all_users": "用户列表",
        "lbl_username": "用户名",
        "lbl_password": "密码",
        "lbl_role": "角色",
        "btn_create_user": "创建"
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
        "qmode_metabo": "⚗️ Metabolite Advanced Analysis",
        "qmode_spatial": "🗺️ Spatial Region Analysis",

        "tab_user_mgmt": "👥 User Mgmt",
        "tab_audit_logs": "📝 Audit Logs",
        "tab_data_crud": "🧬 Data Maintenance",
        "tab_backup": "💾 Backup & Restore",

        "crud_exp_anno": "📝 Expert Annotations (MySQL)",
        "crud_header_core": "🛠️ Core Omics Maintenance (SQLite)",
        "crud_select_db": "Select Database",
        "crud_mode_anno_add": "➕ Add Annotation",
        "crud_mode_anno_manage": "🖊️ Manage Annotations (Edit/Del)",
        "crud_op_create": "➕ Create",
        "crud_op_update": "📝 Update/Extend",
        "crud_op_delete": "🗑️ Delete",
        "crud_op_metabo_map": "🔗 Gene-Metabolite Mapping",
        "crud_op_metabo_exp": "📊 Metabolite Expression",

        "btn_submit": "Submit",
        "btn_save": "Save Changes",
        "btn_delete": "Delete Record",
        "btn_add": "Execute Add",
        "btn_update": "Execute Update",
        "btn_extend": "Extend Column & Update",

        "msg_success": "Operation Successful!",
        "msg_fail": "Operation Failed. Check logs.",
        "msg_confirm_del": "Are you sure you want to delete this?",

        "lbl_target_gene": "Target Gene",
        "lbl_content": "Content",
        "lbl_subtype": "Subtype",
        "lbl_celltype": "CellType",
        "lbl_exp_val": "Value",
        "lbl_sample_id": "Sample ID",
        "lbl_ref_gene": "Reference Gene",
        "lbl_col_gene": "Gene Column",
        "lbl_region": "Region",
        "lbl_metabo": "Metabolite",
        "lbl_kegg": "KEGG Pathway",
        "lbl_note": "Note",
        "lbl_json": "JSON Content",
        "lbl_anno_id": "Annotation ID",
        "lbl_new_sample": "New Sample ID",
        "lbl_auto_col": "⚠️ If gene column missing, table will be altered (ALTER TABLE).",
        "lbl_metabo_mode": "Maintenance Mode",
        "lbl_search_kegg": "Filter by KEGG",
        "lbl_json_example": "Example loaded. Modify based on this:",

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
        "imaging_note": "💡 Note: ROI visualization. JSON defines polygon vertices.",
        "spatial_single_note": "ℹ️ Note: Showing Reference Sample V1 (HER2_Positive).",
        "atac_sim_note": "⚠️ Note: ATAC data shows real + simulated subtypes.",
        "atac_raw_title": "2. Raw Sample Distribution",

        "data_browser": "📚 Data Dictionary",
        "top_genes_list": "🔥 Top Genes (scRNA)",
        "top_metas_list": "🧪 Top Metabolites",
        "input_gene_ph": "Try: FOXA1, ESR1, PKM",
        "input_top_n": "Show Top N",
        "input_sample": "Select Sample",
        "input_metabo": "Select Metabolite",
        "input_region": "Select Region",
        "input_celltype": "Select CellType",
        "header_top_genes": "🔥 Top Expressed Genes",
        "mgmt_create_user": "Create User",
        "mgmt_all_users": "All Users",
        "lbl_username": "Username",
        "lbl_password": "Password",
        "lbl_role": "Role",
        "btn_create_user": "Create"
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

# Imaging Example JSON (For User Convenience)
EXAMPLE_IMAGING_JSON = """{
  "positive": [
    {
      "vertices": [[6000, 14000], [6200, 14000], [6100, 13800]]
    }
  ],
  "negative": []
}"""

if 'logged_in' not in st.session_state:
    st.session_state['logged_in'] = False
    st.session_state['user_role'] = None
    st.session_state['username'] = None
if 'lang' not in st.session_state:
    st.session_state['lang'] = "CN"


def t(key):
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


def get_atac_meta_fallback(sample_id):
    """Fallback hash algorithm if metadata missing"""
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


# --- PDF Generation ---
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

        # === 1. Expert Annotations (Enhanced: Add + Update + Delete) ===
        with st.expander(t('crud_exp_anno'), expanded=True):
            anno_mode = st.radio("Mode", [t('crud_mode_anno_add'), t('crud_mode_anno_manage')], horizontal=True)

            if anno_mode == t('crud_mode_anno_add'):
                c1, c2 = st.columns([1, 2])
                target_gene = c1.text_input(t('lbl_target_gene'), "FOXA1").strip().upper()
                note_content = c2.text_input(t('lbl_content'))
                if st.button(t('btn_submit'), key="btn_anno_add"):
                    conn = auth.get_connection()
                    if conn:
                        cursor = conn.cursor()
                        cursor.execute("INSERT INTO gene_annotations (gene, note, author) VALUES (%s, %s, %s)",
                                       (target_gene, note_content, st.session_state['username']))
                        conn.commit()
                        conn.close()
                        st.success(t('msg_success'))
                        auth.log_action(st.session_state['username'], f"Create Annotation: {target_gene}")
            else:
                # Manage Mode
                conn = auth.get_connection()
                if conn:
                    df_annos = safe_cursor_fetch(conn, "SELECT * FROM gene_annotations ORDER BY created_at DESC")
                    st.dataframe(df_annos, use_container_width=True, height=200)

                    c1, c2, c3 = st.columns([1, 2, 1])
                    selected_id = c1.number_input(t('lbl_anno_id'), min_value=1, step=1)
                    new_note = c2.text_input(t('lbl_content'), key="edit_anno_content")

                    if c3.button(t('btn_save'), key="btn_anno_upd"):
                        cursor = conn.cursor()
                        cursor.execute("UPDATE gene_annotations SET note=%s WHERE id=%s", (new_note, selected_id))
                        if cursor.rowcount > 0:
                            conn.commit()
                            st.success(t('msg_success'))
                        else:
                            st.error(t('msg_fail'))
                        conn.close()

                    if c3.button(t('btn_delete'), key="btn_anno_del"):
                        cursor = conn.cursor()
                        cursor.execute("DELETE FROM gene_annotations WHERE id=%s", (selected_id,))
                        if cursor.rowcount > 0:
                            conn.commit()
                            st.success(t('msg_success'))
                        else:
                            st.error(t('msg_fail'))
                        conn.close()

        st.divider()

        # === 2. Omics Data Maintenance ===
        st.markdown(f"#### {t('crud_header_core')}")
        db_options = ["scRNA", "ATAC", "Metabo", "Spatial", "Imaging"]
        selected_db = st.selectbox(t('crud_select_db'), db_options)

        with st.container(border=True):
            st.markdown(f"**{selected_db} Maintenance**")

            # --- scRNA CRUD (Unchanged as requested) ---
            if selected_db == "scRNA":
                crud_tab1, crud_tab2, crud_tab3 = st.tabs(
                    [t('crud_op_create'), t('crud_op_update'), t('crud_op_delete')])
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
                        except Exception as e:
                            st.error(f"{e}")
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

            # --- ATAC CRUD (Improved: Subtype metadata) ---
            elif selected_db == "ATAC":
                crud_tab1, crud_tab2, crud_tab3 = st.tabs(
                    [t('crud_op_create'), t('crud_op_update'), t('crud_op_delete')])
                conn = sqlite3.connect(DB_PATHS['ATAC'])
                # Ensure metadata table
                conn.execute("CREATE TABLE IF NOT EXISTS sample_metadata (sample TEXT PRIMARY KEY, subtype TEXT)")

                with crud_tab1:  # Create
                    c1, c2, c3 = st.columns(3)
                    new_samp = c1.text_input(t('lbl_new_sample'), key="at_c_s")
                    new_sub = c2.selectbox(t('lbl_subtype'), ["TNBC", "Luminal_A", "HER2_Positive", "Normal"],
                                           key="at_c_sub")
                    ref_gene = c3.text_input(t('lbl_ref_gene'), "FOXA1", key="at_c_g")
                    val = st.number_input(t('lbl_exp_val'), 0.0, key="at_c_v")
                    if st.button(t('btn_add'), key="atac_add"):
                        try:
                            cur = conn.cursor()
                            cur.execute(f"INSERT INTO sample_gene_matrix (sample, {ref_gene}) VALUES (?, ?)",
                                        (new_samp, val))
                            cur.execute("INSERT OR REPLACE INTO sample_metadata (sample, subtype) VALUES (?,?)",
                                        (new_samp, new_sub))
                            conn.commit()
                            st.success(t('msg_success'))
                            auth.log_action(st.session_state['username'], f"Create ATAC Sample: {new_samp}")
                        except Exception as e:
                            st.error(f"{t('msg_fail')}: {e}")

                with crud_tab2:  # Update Value / Add Column
                    st.info(t('lbl_auto_col'))
                    c1, c2, c3 = st.columns(3)
                    samples = pd.read_sql("SELECT sample FROM sample_gene_matrix", conn)['sample'].tolist()
                    tgt_sample = c1.selectbox(t('lbl_sample_id'), samples, key="at_u_s")
                    tgt_gene = c2.text_input(t('lbl_col_gene'), "FOXA1", key="at_u_g").upper()
                    new_val = c3.number_input(t('lbl_new_val'), 0.0, key="at_u_v")

                    if st.button(t('btn_extend'), key="atac_upd"):
                        try:
                            cursor = conn.cursor()
                            try:
                                cursor.execute(f"SELECT {tgt_gene} FROM sample_gene_matrix LIMIT 1")
                            except sqlite3.OperationalError:
                                cursor.execute(f"ALTER TABLE sample_gene_matrix ADD COLUMN {tgt_gene} REAL DEFAULT 0")
                                st.toast(f"Schema Updated: Added column {tgt_gene}")

                            cursor.execute(f"UPDATE sample_gene_matrix SET {tgt_gene} = ? WHERE sample = ?",
                                           (new_val, tgt_sample))
                            conn.commit()
                            st.success(t('msg_success'))
                            auth.log_action(st.session_state['username'], f"Update/Extend ATAC: {tgt_gene}")
                        except Exception as e:
                            st.error(f"{e}")

                with crud_tab3:  # Delete Sample
                    del_sample = st.selectbox(t('lbl_sample_id'), samples, key="atac_del")
                    if st.button(t('btn_delete'), key="atac_del_btn"):
                        conn.execute("DELETE FROM sample_gene_matrix WHERE sample = ?", (del_sample,))
                        conn.execute("DELETE FROM sample_metadata WHERE sample = ?", (del_sample,))
                        conn.commit()
                        st.success(t('msg_success'))
                conn.close()

            # --- Metabo CRUD (Improved: Split into Mapping / Expression) ---
            elif selected_db == "Metabo":
                metabo_mode = st.radio(t('lbl_metabo_mode'), [t('crud_op_metabo_map'), t('crud_op_metabo_exp')],
                                       horizontal=True)
                conn = sqlite3.connect(DB_PATHS['Metabo'])

                if metabo_mode == t('crud_op_metabo_map'):
                    # MAPPING CRUD (Gene <-> Metabolite)
                    crud_tab1, crud_tab2, crud_tab3 = st.tabs(
                        [t('crud_op_create'), t('crud_op_update'), t('crud_op_delete')])
                    with crud_tab1:  # Create Link
                        c1, c2, c3, c4 = st.columns(4)
                        new_g = c1.text_input(t('lbl_target_gene'), key="mt_m_g").upper()
                        new_m = c2.text_input(t('lbl_metabo'), key="mt_m_m")
                        new_k = c3.text_input(t('lbl_kegg'), key="mt_m_k")
                        new_n = c4.text_input(t('lbl_note'), key="mt_m_n")
                        if st.button(t('btn_add'), key="mt_map_add"):
                            conn.execute(
                                "INSERT INTO Gene_Metabolite_Map (Gene, Metabolite, KEGG, Note) VALUES (?,?,?,?)",
                                (new_g, new_m, new_k, new_n))
                            conn.commit()
                            st.success(t('msg_success'))
                    with crud_tab2:  # Update Link Note
                        c1, c2, c3 = st.columns(3)
                        u_g = c1.text_input(t('lbl_target_gene'), key="mt_mu_g").upper()
                        u_m = c2.text_input(t('lbl_metabo'), key="mt_mu_m")
                        u_n = c3.text_input(t('lbl_new_val') + " (Note)", key="mt_mu_n")
                        if st.button(t('btn_update'), key="mt_map_upd"):
                            conn.execute("UPDATE Gene_Metabolite_Map SET Note=? WHERE Gene=? AND Metabolite=?",
                                         (u_n, u_g, u_m))
                            conn.commit()
                            st.success(t('msg_success'))
                    with crud_tab3:  # Delete Link
                        d_g = st.text_input(t('lbl_target_gene'), key="mt_md_g").upper()
                        if st.button(t('btn_delete'), key="mt_map_del"):
                            conn.execute("DELETE FROM Gene_Metabolite_Map WHERE Gene=?", (d_g,))
                            conn.commit()
                            st.success(t('msg_success'))

                else:
                    # EXPRESSION CRUD (Original)
                    crud_tab1, crud_tab2, crud_tab3 = st.tabs(
                        [t('crud_op_create'), t('crud_op_update'), t('crud_op_delete')])
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
                crud_tab1, crud_tab2, crud_tab3 = st.tabs(
                    [t('crud_op_create'), t('crud_op_update'), t('crud_op_delete')])
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

            # --- Imaging CRUD (Improved: Example JSON) ---
            # --- Imaging CRUD (最终修复版: 自动适配表结构) ---
                # --- Imaging CRUD (修复版: 适配 UNIQUE 约束) ---
            elif selected_db == "Imaging":
                conn = sqlite3.connect(DB_PATHS['Imaging'])

                # 尝试建表 (适配组员设计，假设 image_id 存在且唯一)
                # 注意：如果表已存在，这句话会被忽略，不会覆盖已有结构
                conn.execute("""
                                CREATE TABLE IF NOT EXISTS annotations (
                                    id INTEGER PRIMARY KEY AUTOINCREMENT, 
                                    annotation TEXT,
                                    image_id INTEGER UNIQUE
                                )
                            """)

                # 自动补全列 (防止旧文件缺少列)
                try:
                    conn.execute("SELECT image_id FROM annotations LIMIT 1")
                except sqlite3.OperationalError:
                    # 如果没有 image_id，加一列，注意不能设为 UNIQUE (SQLite限制)，只能先加上
                    conn.execute("ALTER TABLE annotations ADD COLUMN image_id INTEGER DEFAULT 0")
                    conn.commit()

                crud_tab1, crud_tab2, crud_tab3 = st.tabs(
                    [t('crud_op_create'), t('crud_op_update'), t('crud_op_delete')])

                with crud_tab1:  # Create
                    st.info(t('lbl_json_example'))

                    c1, c2 = st.columns([1, 3])
                    # 【核心修改】增加输入框，让管理员指定不重复的 ID
                    new_img_id = c1.number_input("Image ID (Unique)", min_value=1, step=1, value=1)

                    default_json = """{
              "positive": [
                {
                  "vertices": [[6000, 14000], [6200, 14000], [6100, 13800]]
                }
              ],
              "negative": []
            }"""
                    json_str = c2.text_area(t('lbl_json'), default_json, height=200)

                    if st.button(t('btn_add'), key="img_add"):
                        try:
                            valid_json = json.loads(json_str)
                            # 使用用户输入的 ID 插入
                            conn.execute("INSERT INTO annotations (annotation, image_id) VALUES (?, ?)",
                                         (json.dumps(valid_json), new_img_id))
                            conn.commit()
                            st.success(t('msg_success'))
                            auth.log_action(st.session_state['username'], f"Create Imaging ROI (ImgID: {new_img_id})")
                        except json.JSONDecodeError:
                            st.error("❌ 格式错误: JSON 无效")
                        except sqlite3.IntegrityError:
                            st.error(f"❌ 提交失败: Image ID {new_img_id} 已存在，请换一个数字！")
                        except Exception as e:
                            st.error(f"❌ 数据库错误: {e}")

                with crud_tab2:  # Update
                    cursor = conn.cursor()
                    try:
                        # 1. 获取所有 ID (刷新列表)
                        df_ids = pd.read_sql("SELECT id, image_id FROM annotations ORDER BY id DESC", conn)

                        if not df_ids.empty:
                            # 2. 构建下拉菜单选项
                            id_map = {}
                            for _, row in df_ids.iterrows():
                                # 确保 ID 是整数
                                rid = int(row['id'])
                                rimg_id = row['image_id']
                                img_id_disp = int(rimg_id) if pd.notna(rimg_id) else "NULL"
                                label = f"ID:{rid} (ImgID:{img_id_disp})"
                                id_map[label] = rid

                            selected_label = st.selectbox("选择要修改的记录", list(id_map.keys()), key="img_upd_sel")

                            # 【核心修复】强制转换为 Python int 类型，解决 numpy.int64 导致的查询失败
                            selected_id = int(id_map[selected_label])

                            # 3. 读取详细数据
                            cursor.execute("SELECT annotation, image_id FROM annotations WHERE id=?", (selected_id,))
                            res = cursor.fetchone()

                            if res:
                                old_json, old_img_id = res
                                safe_img_id = int(old_img_id) if (old_img_id is not None) else 0
                                safe_json = old_json if old_json else "{}"

                                st.markdown("---")
                                c_upd_1, c_upd_2 = st.columns([1, 3])

                                # Image ID 输入框
                                new_img_id_upd = c_upd_1.number_input(
                                    "Image ID (Unique)",
                                    min_value=0, step=1,
                                    value=safe_img_id,
                                    key="img_upd_val_id"
                                )

                                # JSON 输入框
                                new_json = c_upd_2.text_area(
                                    "JSON 数据",
                                    value=safe_json,
                                    height=250,
                                    key="img_upd_json_area"
                                )

                                if st.button("💾 保存修改", key="img_upd_btn"):
                                    try:
                                        json.loads(new_json)  # 校验
                                        conn.execute(
                                            "UPDATE annotations SET annotation=?, image_id=? WHERE id=?",
                                            (new_json, new_img_id_upd, selected_id)
                                        )
                                        conn.commit()
                                        st.success("✅ 修改已保存！")
                                        auth.log_action(st.session_state['username'],
                                                        f"Update Imaging ID: {selected_id}")
                                        time.sleep(0.5)
                                        st.rerun()
                                    except sqlite3.IntegrityError:
                                        st.error(f"❌ 保存失败: Image ID {new_img_id_upd} 已被占用，请换一个数字。")
                                    except json.JSONDecodeError:
                                        st.error("❌ JSON 格式错误")
                            else:
                                # 增加详细调试信息，万一再报错知道原因
                                st.error(f"❌ 未找到 ID={selected_id} 的记录。建议刷新页面重试。")
                        else:
                            st.info("📭 数据库目前是空的，请先去 [新增] 标签页添加数据。")
                    except Exception as e:
                        st.error(f"❌ 系统错误: {e}")

                with crud_tab3:  # Delete
                    df_ids = pd.read_sql("SELECT id, image_id FROM annotations", conn)
                    if not df_ids.empty:
                        id_map = {f"ID:{row['id']} (ImgID:{row['image_id']})": row['id'] for _, row in
                                  df_ids.iterrows()}
                        del_label = st.selectbox("选择要删除的记录", list(id_map.keys()), key="img_del_sel")
                        id_to_del = id_map[del_label]

                        if st.button(t('btn_delete'), key="img_del"):
                            conn.execute("DELETE FROM annotations WHERE id=?", (id_to_del,))
                            conn.commit()
                            st.success(t('msg_success'))
                            auth.log_action(st.session_state['username'], f"Delete Imaging ID: {id_to_del}")
                            st.rerun()
                    else:
                        st.info("无数据可删除")

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
                # Optimized: Try fetching real metadata first
                meta_dict = {}
                conn_atac = sqlite3.connect(DB_PATHS['ATAC'])
                try:
                    meta_df = pd.read_sql("SELECT sample, subtype FROM sample_metadata", conn_atac)
                    meta_dict = dict(zip(meta_df['sample'], meta_df['subtype']))
                except:
                    pass
                conn_atac.close()

                # Apply metadata
                df['Subtype'] = df['sample'].apply(lambda x: meta_dict.get(x, get_atac_meta_fallback(x)))

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

    # --- Advanced Modes ---
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
                st.bar_chart(df.head(top_n), x="Gene", y="MeanExpr")
                st.dataframe(df, use_container_width=True)

    elif mode_key == "atac_advanced":
        st.subheader(t('qmode_atac'))
        samples = get_distinct_values("ATAC", "sample_gene_matrix", "sample")
        c1, c2 = st.columns(2)
        sid = c1.selectbox(t('input_sample'), samples)
        top_n = c2.slider(t('input_top_n'), 3, 100, 5)
        if sid:
            df = run_sqlite_query("ATAC", f"SELECT * FROM sample_gene_matrix WHERE sample = '{sid}'")
            if df is not None:
                df_t = df.drop(columns=['sample']).T
                df_t.columns = ['Openness']
                df_t = df_t.sort_values(by='Openness', ascending=False).head(top_n)
                st.bar_chart(df_t)

    elif mode_key == "metabo_advanced":
        st.subheader(t('qmode_metabo'))
        # Enhanced Metabo Query
        c1, c2 = st.columns(2)
        metas = get_distinct_values("Metabo", "Gene_Metabolite_Map", "Metabolite")
        keggs = get_distinct_values("Metabo", "Gene_Metabolite_Map", "KEGG")

        filter_type = c1.radio("Filter By", ["Metabolite", "KEGG"], horizontal=True)

        if filter_type == "Metabolite":
            m = c2.selectbox(t('input_metabo'), metas)
            if m:
                m_safe = m.replace("'", "''")
                df = run_sqlite_query("Metabo",
                                      f"SELECT Subtype, Expression_Level FROM Metabolite_Expression WHERE Metabolite = '{m_safe}'")
                st.markdown(f"#### Expression of {m}")
                if df is not None: st.bar_chart(df, x="Subtype", y="Expression_Level")
        else:
            k = c2.selectbox(t('lbl_search_kegg'), keggs)
            if k:
                k_safe = k.replace("'", "''")
                df_genes = run_sqlite_query("Metabo",
                                            f"SELECT Gene, Metabolite, Note FROM Gene_Metabolite_Map WHERE KEGG = '{k_safe}'")
                st.markdown(f"#### Pathway: {k}")
                st.dataframe(df_genes, use_container_width=True)

    elif mode_key == "spatial_advanced":
        st.subheader(t('qmode_spatial'))
        regions = get_distinct_values("Spatial", "Table_SpatialExpression", "Region")
        c1, c2 = st.columns(2)
        r = c1.selectbox(t('input_region'), regions)
        top_n = c2.slider(t('input_top_n'), 3, 100, 5)
        if r:
            df = run_sqlite_query("Spatial",
                                  f"SELECT Gene, Avg_Expression FROM Table_SpatialExpression WHERE Region = '{r}' ORDER BY Avg_Expression DESC LIMIT {top_n}")
            if df is not None: st.bar_chart(df, x="Gene", y="Avg_Expression")


if __name__ == "__main__":
    if st.session_state['logged_in']:
        query_ui()
    else:
        login_ui()