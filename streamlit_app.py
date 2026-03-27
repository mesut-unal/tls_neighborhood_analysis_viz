import os
import glob
import io
from collections import OrderedDict
from typing import List, Dict, Tuple

import streamlit as st
from pdf2image import convert_from_path


def sanitize_for_filename(name: str) -> str:
    out = name.replace(" ", "_")
    for ch in [":", "/", "\\"]:
        out = out.replace(ch, "")
    return out


def folder_templates(postfix: str) -> Dict[str, str]:
    return {
        "pval_barplots": f"tls_pval_barplots_{postfix}",
        "gmm_components": f"tail_3sig_byMedian_median_plus_3mad_{postfix}",
        "violin": f"tls_violin_by_tls_cluster_{postfix}",
        "radius_profiles": f"tls_radius_profiles_grouped_variation_gmmMedian_allPatients_{postfix}",
    }


def resolve_path(base_dir: str, *parts: str) -> str:
    return os.path.join(base_dir, *parts)


def file_exists(path: str) -> bool:
    return os.path.isfile(path)


def show_pdf_as_image(container, path: str, caption: str):
    try:
        pages = convert_from_path(path, dpi=150, first_page=1, last_page=1)
        buf = io.BytesIO()
        pages[0].save(buf, format="PNG")
        buf.seek(0)

        container.image(buf.getvalue(), caption=caption, width="stretch")

        with open(path, "rb") as f:
            pdf_bytes = f.read()
        container.download_button(
            label=f"⬇️ Download PDF",
            data=pdf_bytes,
            file_name=os.path.basename(path),
            mime="application/pdf",
            key=f"dl_{path}",
        )
    except Exception as e:
        container.error(f"Could not render {os.path.basename(path)}: {e}")


def find_markers_from_patterns(patterns: List[str], base_dir: str, gmm_dir: str) -> List[str]:
    found = set()
    for pat in patterns:
        for grp in ("CLR", "DII"):
            suffix = f"_{grp}_gmm_components.pdf"
            candidates = [
                resolve_path(base_dir, gmm_dir, f"{pat}{suffix}"),
                resolve_path(base_dir, gmm_dir, f"{sanitize_for_filename(pat)}{suffix}"),
            ]
            for c in candidates:
                for path in glob.glob(c):
                    fname = os.path.basename(path)
                    end = f"_{grp}_gmm_components.pdf"
                    if fname.endswith(end):
                        base = fname[: -len(end)]
                        found.add(base)
                        break
    return sorted(found)


def image_paths_for_marker(marker: str, base_dir: str, dirs: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    m_raw = marker
    m_san = sanitize_for_filename(marker)

    def try_paths(folder: str, prefix: str, suffix: str) -> Tuple[str, str]:
        p_clr = resolve_path(base_dir, folder, f"{prefix}_CLR_{suffix}")
        p_dii = resolve_path(base_dir, folder, f"{prefix}_DII_{suffix}")
        if not file_exists(p_clr):
            p_clr_alt = resolve_path(base_dir, folder, f"{m_san}_CLR_{suffix}")
            if file_exists(p_clr_alt):
                p_clr = p_clr_alt
        if not file_exists(p_dii):
            p_dii_alt = resolve_path(base_dir, folder, f"{m_san}_DII_{suffix}")
            if file_exists(p_dii_alt):
                p_dii = p_dii_alt
        return p_clr, p_dii

    gmm_clr, gmm_dii = try_paths(dirs["gmm_components"], m_raw, "gmm_components.pdf")
    vio_clr, vio_dii = try_paths(dirs["violin"], m_raw, "violin.pdf")
    rad_clr, rad_dii = try_paths(dirs["radius_profiles"], m_raw, "radius_profile.pdf")

    return {
        "gmm_components": {"CLR": gmm_clr, "DII": gmm_dii},
        "violin": {"CLR": vio_clr, "DII": vio_dii},
        "radius": {"CLR": rad_clr, "DII": rad_dii},
    }


st.set_page_config(page_title="TLS Results Viewer", layout="wide")
st.title("TLS Results Viewer")

base_dir = "."
postfix = "26MAR2026"

MARKER_GROUPS: "OrderedDict[str, List[str]]" = OrderedDict({
    "EGFR marker": ["EGFR*"],
    "Macrophages markers": ["CD11b*", "CD11c*", "CD68*", "CD163*"],
    "NK markers": ["CD56*", "CD57*"],
    "Phagocytes markers": ["BCL-2*"],
    "Proliferation markers": ["Ki67*"],
    "Stroma markers": ["aSMA*", "beta-catenin*", "MMP9*", "MMP12*", "CD44*", "GFAP*", "CD15*", "Collagen*"],
    "T cell markers": ["CD194*", "GATA3*", "PD1*", "ICOS*", "LAG3*", "FOXP3*"],
    "Vasculatures markers": ["CD31*", "CD34*", "CD44*"],
})

DIRS = folder_templates(postfix)

# ── p-value barplots ──────────────────────────────────────────────────────────
st.header("Inside vs. outside TLS structure marker intensity differences")

pval_dir = resolve_path(base_dir, DIRS["pval_barplots"])

clr_files = [
    "CLR_nascent_pvalue_barplot_kruskal.pdf",
    "CLR_intermediate_pvalue_barplot_kruskal.pdf",
    "CLR_mature_pvalue_barplot_kruskal.pdf",
]
dii_files = [
    "DII_nascent_pvalue_barplot_kruskal.pdf",
    "DII_intermediate_pvalue_barplot_kruskal.pdf",
]

cols = st.columns(3)
for i, fname in enumerate(clr_files):
    path = resolve_path(pval_dir, fname)
    if file_exists(path):
        show_pdf_as_image(cols[i], path, fname)
    else:
        cols[i].info(f"Missing: {path}")

cols = st.columns(3)
for i, fname in enumerate(dii_files):
    path = resolve_path(pval_dir, fname)
    if file_exists(path):
        show_pdf_as_image(cols[i], path, fname)
    else:
        cols[i].info(f"Missing: {path}")

st.divider()

# ── per-marker panels ─────────────────────────────────────────────────────────
st.header("Marker categories")

for category, patterns in MARKER_GROUPS.items():
    with st.expander(label=category, expanded=False):
        st.markdown(f"<h3 style='font-size:20px; font-weight:bold; margin-top:0'>{category}</h3>", unsafe_allow_html=True)
        markers = find_markers_from_patterns(patterns, base_dir, DIRS["gmm_components"])
        if not markers:
            st.warning("No markers found for patterns: " + ", ".join(patterns))
            continue

        for marker in markers:
            st.subheader(marker)
            paths = image_paths_for_marker(marker, base_dir, DIRS)
            col1, col2, col3 = st.columns(3, gap="large")

            with col1:
                st.caption("Prob. density vs. expression level")
                for grp in ("CLR", "DII"):
                    p = paths["gmm_components"][grp]
                    if file_exists(p):
                        show_pdf_as_image(st, p, f"{marker} {grp} gmm_components")

            with col2:
                st.caption("Marker intensity across maturation levels")
                for grp in ("CLR", "DII"):
                    p = paths["violin"][grp]
                    if file_exists(p):
                        show_pdf_as_image(st, p, f"{marker} {grp} violin")

            with col3:
                st.caption("Average intensity at different radii")
                for grp in ("CLR", "DII"):
                    p = paths["radius"][grp]
                    if file_exists(p):
                        show_pdf_as_image(st, p, f"{marker} {grp} radius_profile")

            st.markdown("---")

# ── Marker vs. Cell Types ─────────────────────────────────────────────────────
st.header("Marker vs. Cell Types")
celltype_folder = resolve_path(base_dir, f"marker_vs_celltype_two_panel_{postfix}", "cell_type")

if os.path.isdir(celltype_folder):
    imgs = sorted(glob.glob(os.path.join(celltype_folder, "*.pdf")))
    if imgs:
        for i in range(0, len(imgs), 3):
            row = imgs[i:i+3]
            cols = st.columns(len(row))
            for c, p in zip(cols, row):
                show_pdf_as_image(c, p, os.path.basename(p))
    else:
        st.info(f"No .pdf files found in {celltype_folder}")
else:
    st.info(f"Folder not found: {celltype_folder}")

# ── Marker vs. Cell States ────────────────────────────────────────────────────
st.header("Marker vs. Cell States")
cellstate_folder = resolve_path(base_dir, f"marker_vs_celltype_two_panel_{postfix}", "cell_state")

if os.path.isdir(cellstate_folder):
    imgs = sorted(glob.glob(os.path.join(cellstate_folder, "*.pdf")))
    if imgs:
        for i in range(0, len(imgs), 3):
            row = imgs[i:i+3]
            cols = st.columns(len(row))
            for c, p in zip(cols, row):
                show_pdf_as_image(c, p, os.path.basename(p))
    else:
        st.info(f"No .pdf files found in {cellstate_folder}")
else:
    st.info(f"Folder not found: {cellstate_folder}")