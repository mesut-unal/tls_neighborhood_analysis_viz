import os
import glob
from collections import OrderedDict
from typing import List, Dict, Tuple

import streamlit as st

def sanitize_for_filename(name: str) -> str:
    out = name.replace(" ", "_")
    for ch in [":", "/", "\\"]:
        out = out.replace(ch, "")
    return out

def folder_templates(postfix: str) -> Dict[str, str]:
    return {
        "pval_barplots": f"tls_pval_barplots_{postfix}",
        "gmm_components": f"tail_3sig_byMedian_{postfix}_median_plus_3mad",
        "violin": f"tls_violin_by_tls_cluster_{postfix}",
        "radius_profiles": f"tls_radius_profiles_grouped_variation_gmmMedian_allPatients_{postfix}_inwardLimit2",
    }

def resolve_path(base_dir: str, *parts: str) -> str:
    return os.path.join(base_dir, *parts)

def file_exists(path: str) -> bool:
    return os.path.isfile(path)

def find_markers_from_patterns(patterns: List[str], base_dir: str, gmm_dir: str) -> List[str]:
    found = set()
    for pat in patterns:
        for grp in ("CRL", "DII"):
            suffix = f"_{grp}_gmm_components.png"
            candidates = [
                resolve_path(base_dir, gmm_dir, f"{pat}{suffix}"),
                resolve_path(base_dir, gmm_dir, f"{sanitize_for_filename(pat)}{suffix}"),
            ]
            for c in candidates:
                for path in glob.glob(c):
                    fname = os.path.basename(path)
                    end = f"_{grp}_gmm_components.png"
                    if fname.endswith(end):
                        base = fname[: -len(end)]
                        found.add(base)
                        break
    return sorted(found)

def image_paths_for_marker(marker: str, base_dir: str, dirs: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    m_raw = marker
    m_san = sanitize_for_filename(marker)

    def try_paths(folder: str, prefix: str, suffix: str) -> Tuple[str, str]:
        p_crl = resolve_path(base_dir, folder, f"{prefix}_CRL_{suffix}")
        p_dii = resolve_path(base_dir, folder, f"{prefix}_DII_{suffix}")
        if not file_exists(p_crl):
            p_crl_alt = resolve_path(base_dir, folder, f"{m_san}_CRL_{suffix}")
            if file_exists(p_crl_alt):
                p_crl = p_crl_alt
        if not file_exists(p_dii):
            p_dii_alt = resolve_path(base_dir, folder, f"{m_san}_DII_{suffix}")
            if file_exists(p_dii_alt):
                p_dii = p_dii_alt
        return p_crl, p_dii

    gmm_crl, gmm_dii = try_paths(dirs["gmm_components"], m_raw, "gmm_components.png")
    vio_crl, vio_dii = try_paths(dirs["violin"], m_raw, "violin.png")
    rad_crl, rad_dii = try_paths(dirs["radius_profiles"], m_raw, "radius_profile.png")

    return {
        "gmm_components": {"CRL": gmm_crl, "DII": gmm_dii},
        "violin": {"CRL": vio_crl, "DII": vio_dii},
        "radius": {"CRL": rad_crl, "DII": rad_dii},
    }

st.set_page_config(page_title="TLS Results Viewer", layout="wide")

st.title("TLS Results Viewer")

base_dir = "."
postfix = "aug14"

MARKER_GROUPS: "OrderedDict[str, List[str]]" = OrderedDict({
    "EGFR marker": ["EGFR*"],
    "Macrophages markers": ["CD11b*", "CD11c*"],
    "NK markers": ["CD56*", "CD57*"],
    "Phagocytes markers": ["BCL-2*"],
    "Proliferation markers": ["Ki67*"],
    "Stroma markers": ["aSMA*", "beta-catenin*", "MMP9*", "MMP12*", "CD44*", "GFAP*", "CD15*"],
    "T cell markers": ["CD194*", "GATA3*", "PD1*", "ICOS*", "LAG3*", "FOXP3*"],
    "Vasculatures markers": ["CD31*", "CD34*", "CD44*"],
})

DIRS = folder_templates(postfix)

st.header("Inside vs. outside TLS structure marker intensity differences")

pval_dir = resolve_path(base_dir, DIRS["pval_barplots"])

crl_files = [
    "CRL_nascent_pvalue_barplot_kruskal.png",
    "CRL_intermediate_pvalue_barplot_kruskal.png",
    "CRL_mature_pvalue_barplot_kruskal.png",
]

dii_files = [
    "DII_nascent_pvalue_barplot_kruskal.png",
    "DII_intermediate_pvalue_barplot_kruskal.png",
]

cols = st.columns(3)
for i, fname in enumerate(crl_files):
    path = resolve_path(pval_dir, fname)
    if file_exists(path):
        cols[i].image(path, caption=fname, use_container_width=True)
    else:
        cols[i].info(f"Missing: {path}")

cols = st.columns(3)
for i, fname in enumerate(dii_files):
    path = resolve_path(pval_dir, fname)
    if file_exists(path):
        cols[i].image(path, caption=fname, use_container_width=True)
    else:
        cols[i].info(f"Missing: {path}")

st.divider()

st.header("Marker categories")

for category, patterns in MARKER_GROUPS.items():
    # Using unsafe_allow_html to style the expander header text without showing HTML tags
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
                gmm_crl = paths["gmm_components"]["CRL"]
                gmm_dii = paths["gmm_components"]["DII"]
                if file_exists(gmm_crl):
                    st.image(gmm_crl, caption=f"{marker} CRL gmm_components", use_container_width=True)
                if file_exists(gmm_dii):
                    st.image(gmm_dii, caption=f"{marker} DII gmm_components", use_container_width=True)

            with col2:
                st.caption("Marker intensity accross maturation levels")
                vio_crl = paths["violin"]["CRL"]
                vio_dii = paths["violin"]["DII"]
                if file_exists(vio_crl):
                    st.image(vio_crl, caption=f"{marker} CRL violin", use_container_width=True)
                if file_exists(vio_dii):
                    st.image(vio_dii, caption=f"{marker} DII violin", use_container_width=True)

            with col3:
                st.caption("Average intensity at different radii")
                rad_crl = paths["radius"]["CRL"]
                rad_dii = paths["radius"]["DII"]
                if file_exists(rad_crl):
                    st.image(rad_crl, caption=f"{marker} CRL radius_profile", use_container_width=True)
                if file_exists(rad_dii):
                    st.image(rad_dii, caption=f"{marker} DII radius_profile", use_container_width=True)

            st.markdown("---")
