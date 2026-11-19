// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

use crate::model::{
    ReceptorProfile,
    AtlasScale, BrainRegion, ConnectomeAtlas, Hemisphere, PathwayEdge, PathwayModality,
};
use std::collections::BTreeSet;

pub fn all_presets() -> Vec<ConnectomeAtlas> {
    vec![
        human_mni152_atlas(),
        mouse_ccfv3_atlas(),
        drosophila_flywire_atlas(),
        celegans_wormatlas(),
    ]
}

pub fn load_preset_by_id(id: &str) -> Option<ConnectomeAtlas> {
    all_presets().into_iter().find(|p| p.id == id)
}

fn reg(
    id: &str,
    name: &str,
    network: &str,
    hemisphere: Hemisphere,
    coord: [f32; 3],
    radius_mm: f32,
    color_hex: &str,
    citation_url: &str,
    notes: &str,
) -> BrainRegion {
    BrainRegion {
        id: id.to_string(),
        name: name.to_string(),
        network: network.to_string(),
        hemisphere,
        coord,
        radius_mm,
        color_hex: color_hex.to_string(),
        citation_url: citation_url.to_string(),
        notes: notes.to_string(),
        what_it_does: notes.to_string(),
        inputs_prose: "Upstream afferent pathways from connected regions in this connectome.".to_string(),
        outputs_prose: "Downstream efferent projections to target regions in this connectome.".to_string(),
        receptors: infer_receptor_profile(id, network),
    }
}

fn infer_receptor_profile(id: &str, network: &str) -> ReceptorProfile {
    if id.contains("Caud") || id.contains("Put") || id == "CP_L" || id.contains("MB_Lobes") {
        ReceptorProfile { da: 0.95, serotonin: 0.50, ach: 0.89, glutamate: 0.66, gaba: 0.74 }
    } else if id.contains("Hipp") || id == "CA1_L" || id.contains("Amyg") || id == "BLA_L" {
        ReceptorProfile { da: 0.62, serotonin: 0.92, ach: 0.80, glutamate: 0.86, gaba: 0.65 }
    } else if id == "DD02" || id == "EB" || network.contains("Visual") || id.contains("CBL") {
        ReceptorProfile { da: 0.30, serotonin: 0.44, ach: 0.54, glutamate: 0.68, gaba: 0.94 }
    } else {
        ReceptorProfile { da: 0.58, serotonin: 0.68, ach: 0.72, glutamate: 0.84, gaba: 0.66 }
    }
}

fn edge(
    source: &str,
    target: &str,
    weight: f64,
    tract_name: &str,
    modality: PathwayModality,
    directed: bool,
    synapse_count_est: Option<u64>,
) -> PathwayEdge {
    PathwayEdge {
        source: source.to_string(),
        target: target.to_string(),
        weight,
        tract_name: tract_name.to_string(),
        modality,
        directed,
        synapse_count_est,
    }
}

/// 1. Human Macro-Connectome in MNI152 stereotaxic space (mm).
pub fn human_mni152_atlas() -> ConnectomeAtlas {
    let url_conn = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11939878/";
    let url_rn = "https://pmc.ncbi.nlm.nih.gov/articles/PMC9491841/";
    let url_sn = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8859259/";
    let url_vis = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11922427/";
    let url_lang = "https://www.frontiersin.org/journals/human-neuroscience/articles/10.3389/fnhum.2023.1095413/full";
    let url_amyg = "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0095740";
    let url_hipp = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7375100/";
    let url_bg = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6444302/";

    let regions = vec![
        // Default Mode Network (DMN)
        reg("DMN.MPFC", "Medial Prefrontal Cortex", "Default Mode", Hemisphere::Midline, [1.0, 55.0, -3.0], 5.5, "#38bdf8", url_conn, "Anterior midline DMN core hub in MNI152 space"),
        reg("DMN.PCC", "Posterior Cingulate / Precuneus", "Default Mode", Hemisphere::Midline, [1.0, -61.0, 38.0], 6.0, "#38bdf8", url_conn, "Posterior midline DMN integration hub"),
        reg("DMN.LP_L", "Left Lateral Parietal (TPJ)", "Default Mode", Hemisphere::Left, [-39.0, -77.0, 33.0], 4.8, "#38bdf8", url_conn, "Left inferior parietal / temporoparietal DMN node"),
        reg("DMN.LP_R", "Right Lateral Parietal (TPJ)", "Default Mode", Hemisphere::Right, [47.0, -67.0, 29.0], 4.8, "#38bdf8", url_conn, "Right inferior parietal / temporoparietal DMN node"),
        // Sensorimotor Network (SMN)
        reg("SMN.Lat_L", "Left Lateral Sensorimotor", "Sensorimotor", Hemisphere::Left, [-55.0, -12.0, 29.0], 5.0, "#34d399", url_conn, "Left precentral/postcentral sensorimotor cortex"),
        reg("SMN.Lat_R", "Right Lateral Sensorimotor", "Sensorimotor", Hemisphere::Right, [56.0, -10.0, 29.0], 5.0, "#34d399", url_conn, "Right precentral/postcentral sensorimotor cortex"),
        reg("SMN.Sup", "Superior Sensorimotor (Paracentral)", "Sensorimotor", Hemisphere::Midline, [0.0, -31.0, 67.0], 5.2, "#34d399", url_conn, "Dorsal midline primary motor/somatosensory cortex"),
        // Visual Network (VIS)
        reg("VIS.Med", "Medial Visual (Calcarine)", "Visual", Hemisphere::Midline, [2.0, -79.0, 12.0], 5.2, "#a78bfa", url_vis, "Primary striate/calcarine visual cortex"),
        reg("VIS.Occ", "Occipital Visual Pole", "Visual", Hemisphere::Midline, [0.0, -93.0, -4.0], 4.8, "#a78bfa", url_vis, "Posterior occipital visual pole"),
        reg("VIS.Lat_L", "Left Lateral Visual", "Visual", Hemisphere::Left, [-37.0, -79.0, 10.0], 4.6, "#a78bfa", url_vis, "Left extrastriate lateral visual area"),
        reg("VIS.Lat_R", "Right Lateral Visual", "Visual", Hemisphere::Right, [38.0, -72.0, 13.0], 4.6, "#a78bfa", url_vis, "Right extrastriate lateral visual area"),
        // Salience Network (SAL)
        reg("SAL.ACC", "Dorsal Anterior Cingulate", "Salience", Hemisphere::Midline, [0.0, 22.0, 35.0], 5.2, "#fb923c", url_sn, "Midline dorsal ACC salience & executive switch hub"),
        reg("SAL.AIns_L", "Left Anterior Insula", "Salience", Hemisphere::Left, [-44.0, 13.0, 1.0], 4.6, "#fb923c", url_sn, "Left fronto-insular salience hub"),
        reg("SAL.AIns_R", "Right Anterior Insula", "Salience", Hemisphere::Right, [47.0, 14.0, 0.0], 4.6, "#fb923c", url_sn, "Right fronto-insular salience hub"),
        reg("SAL.RPFC_L", "Left Rostral Prefrontal", "Salience", Hemisphere::Left, [-32.0, 45.0, 27.0], 4.4, "#fb923c", url_rn, "Left anterior prefrontal salience parcel"),
        reg("SAL.RPFC_R", "Right Rostral Prefrontal", "Salience", Hemisphere::Right, [32.0, 46.0, 27.0], 4.4, "#fb923c", url_rn, "Right anterior prefrontal salience parcel"),
        reg("SAL.SMG_L", "Left Supramarginal Gyrus", "Salience", Hemisphere::Left, [-60.0, -39.0, 31.0], 4.4, "#fb923c", url_rn, "Left ventral parietal salience node"),
        reg("SAL.SMG_R", "Right Supramarginal Gyrus", "Salience", Hemisphere::Right, [62.0, -35.0, 32.0], 4.4, "#fb923c", url_rn, "Right ventral parietal salience node"),
        // Dorsal Attention Network (DAN)
        reg("DAN.FEF_L", "Left Frontal Eye Field", "Dorsal Attention", Hemisphere::Left, [-27.0, -9.0, 64.0], 4.5, "#2dd4bf", url_rn, "Left dorsal premotor attentional eye field"),
        reg("DAN.FEF_R", "Right Frontal Eye Field", "Dorsal Attention", Hemisphere::Right, [30.0, -6.0, 64.0], 4.5, "#2dd4bf", url_rn, "Right dorsal premotor attentional eye field"),
        reg("DAN.IPS_L", "Left Intraparietal Sulcus", "Dorsal Attention", Hemisphere::Left, [-39.0, -43.0, 52.0], 4.6, "#2dd4bf", url_rn, "Left superior parietal visuospatial attention node"),
        reg("DAN.IPS_R", "Right Intraparietal Sulcus", "Dorsal Attention", Hemisphere::Right, [39.0, -42.0, 54.0], 4.6, "#2dd4bf", url_rn, "Right superior parietal visuospatial attention node"),
        // Frontoparietal Control Network (FPN)
        reg("FPN.LPFC_L", "Left Dorsolateral Prefrontal", "Frontoparietal", Hemisphere::Left, [-43.0, 33.0, 28.0], 4.9, "#facc15", url_rn, "Left lateral prefrontal cognitive control hub"),
        reg("FPN.LPFC_R", "Right Dorsolateral Prefrontal", "Frontoparietal", Hemisphere::Right, [41.0, 38.0, 30.0], 4.9, "#facc15", url_conn, "Right lateral prefrontal cognitive control hub"),
        reg("FPN.PPC_L", "Left Posterior Parietal", "Frontoparietal", Hemisphere::Left, [-46.0, -58.0, 49.0], 4.7, "#facc15", url_conn, "Left posterior parietal executive hub"),
        reg("FPN.PPC_R", "Right Posterior Parietal", "Frontoparietal", Hemisphere::Right, [52.0, -52.0, 45.0], 4.7, "#facc15", url_conn, "Right posterior parietal executive hub"),
        // Perisylvian Language Network (LANG)
        reg("LANG.Broca_L", "Left IFG (Broca's Area)", "Language", Hemisphere::Left, [-51.0, 27.0, 18.0], 4.7, "#f472b6", url_lang, "Left inferior frontal speech production hub"),
        reg("LANG.Wernicke_L", "Left pSTG (Wernicke's Area)", "Language", Hemisphere::Left, [-51.0, -51.0, 30.0], 4.7, "#f472b6", url_lang, "Left posterior superior temporal comprehension hub"),
        // Subcortical & Limbic Nuclei (SUBC)
        reg("SUBC.Amyg_L", "Left Amygdala", "Subcortical & Limbic", Hemisphere::Left, [-20.0, -4.0, -15.0], 4.0, "#f87171", url_amyg, "Left medial temporal amygdalar complex"),
        reg("SUBC.Amyg_R", "Right Amygdala", "Subcortical & Limbic", Hemisphere::Right, [22.0, -2.0, -15.0], 4.0, "#f87171", url_amyg, "Right medial temporal amygdalar complex"),
        reg("SUBC.Hipp_L", "Left Hippocampus", "Subcortical & Limbic", Hemisphere::Left, [-36.0, -9.0, -24.0], 4.4, "#f87171", url_hipp, "Left hippocampal episodic memory formation"),
        reg("SUBC.Thal_L", "Left Thalamus", "Subcortical & Limbic", Hemisphere::Left, [-16.0, -23.0, -4.0], 4.8, "#f87171", url_bg, "Left diencephalic thalamocortical relay nucleus"),
        reg("SUBC.Caud_L", "Left Caudate Nucleus", "Subcortical & Limbic", Hemisphere::Left, [-16.0, 4.0, 14.0], 4.3, "#f87171", url_bg, "Left dorsal striatal associative nucleus"),
        reg("SUBC.Put_L", "Left Putamen", "Subcortical & Limbic", Hemisphere::Left, [-28.0, -2.0, -6.0], 4.5, "#f87171", "https://pmc.ncbi.nlm.nih.gov/articles/PMC6869828/", "Left dorsal striatal motor execution nucleus receiving sensorimotor corticostriatal inputs and projecting to globus pallidus"),
        reg("SUBC.Pall_L", "Left Globus Pallidus", "Subcortical & Limbic", Hemisphere::Left, [-20.0, -6.0, -1.0], 4.0, "#f87171", url_bg, "Left basal ganglia inhibitory output brake gating thalamocortical motor & executive loops"),
        reg("SUBC.RN_L", "Left Red Nucleus (Midbrain)", "Subcortical & Limbic", Hemisphere::Left, [-8.0, -19.0, -8.0], 3.7, "#f87171", url_rn, "Left midbrain tegmental relay integrating corticorubral and crossed cerebellar dentate signals"),
        reg("BSTEM.Pons", "Pontine Nuclei (Brainstem)", "Cerebellar & Brainstem", Hemisphere::Midline, [-8.0, -42.0, -28.0], 4.6, "#38bdf8", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12635369/", "Ventral brainstem pontine relay translating prefrontal and motor cortical plans into cerebellar mossy fibers"),
        reg("CBL.LobVI_R", "Right Cerebellum Lobule VI", "Cerebellar & Brainstem", Hemisphere::Right, [20.0, -50.0, -22.0], 4.6, "#60a5fa", "https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2021.778201/full", "Superior paravermal cerebellar lobule computing millisecond sensorimotor forward models"),
        reg("CBL.Crus1_L", "Left Cerebellum Crus I (HVIIA)", "Cerebellar & Brainstem", Hemisphere::Left, [-44.0, -42.0, -34.0], 5.0, "#60a5fa", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12635369/", "Left lateral neocerebellar cognitive hemisphere coupled with executive and default-mode networks"),
        reg("CBL.Crus1_R", "Right Cerebellum Crus I (HVIIA)", "Cerebellar & Brainstem", Hemisphere::Right, [46.0, -58.0, -34.0], 5.0, "#60a5fa", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12635369/", "Right lateral neocerebellar cognitive hemisphere coupled with left prefrontal and language networks"),
    ];

    let edges = vec![
        // Cingulum Bundle
        edge("DMN.MPFC", "DMN.PCC", 0.92, "Cingulum Bundle", PathwayModality::Association, false, None),
        edge("DMN.MPFC", "SAL.ACC", 0.84, "Cingulum Bundle", PathwayModality::Association, false, None),
        edge("SAL.ACC", "DMN.PCC", 0.86, "Cingulum Bundle", PathwayModality::Association, false, None),
        edge("DMN.PCC", "SUBC.Hipp_L", 0.78, "Cingulum Bundle", PathwayModality::Association, false, None),
        edge("DMN.PCC", "DMN.LP_L", 0.81, "Cingulum Bundle", PathwayModality::Association, false, None),
        edge("DMN.PCC", "DMN.LP_R", 0.80, "Cingulum Bundle", PathwayModality::Association, false, None),
        // Corpus Callosum (Interhemispheric Commissural)
        edge("FPN.LPFC_L", "FPN.LPFC_R", 0.88, "Corpus Callosum", PathwayModality::Commissural, false, None),
        edge("SAL.RPFC_L", "SAL.RPFC_R", 0.83, "Corpus Callosum", PathwayModality::Commissural, false, None),
        edge("DAN.FEF_L", "DAN.FEF_R", 0.85, "Corpus Callosum", PathwayModality::Commissural, false, None),
        edge("SMN.Lat_L", "SMN.Lat_R", 0.87, "Corpus Callosum", PathwayModality::Commissural, false, None),
        edge("DAN.IPS_L", "DAN.IPS_R", 0.82, "Corpus Callosum", PathwayModality::Commissural, false, None),
        edge("FPN.PPC_L", "FPN.PPC_R", 0.84, "Corpus Callosum", PathwayModality::Commissural, false, None),
        edge("VIS.Lat_L", "VIS.Lat_R", 0.89, "Corpus Callosum", PathwayModality::Commissural, false, None),
        edge("DMN.LP_L", "DMN.LP_R", 0.76, "Corpus Callosum", PathwayModality::Commissural, false, None),
        edge("SAL.AIns_L", "SAL.AIns_R", 0.74, "Corpus Callosum", PathwayModality::Commissural, false, None),
        // Arcuate Fasciculus (Language Loop)
        edge("LANG.Broca_L", "LANG.Wernicke_L", 0.91, "Arcuate Fasciculus", PathwayModality::Association, false, None),
        edge("LANG.Broca_L", "SAL.SMG_L", 0.77, "Arcuate Fasciculus", PathwayModality::Association, false, None),
        edge("SAL.SMG_L", "LANG.Wernicke_L", 0.79, "Arcuate Fasciculus", PathwayModality::Association, false, None),
        edge("LANG.Broca_L", "FPN.LPFC_L", 0.73, "Arcuate Fasciculus", PathwayModality::Association, false, None),
        // Superior Longitudinal Fasciculus (SLF)
        edge("FPN.LPFC_L", "FPN.PPC_L", 0.89, "Superior Longitudinal Fasciculus", PathwayModality::Association, false, None),
        edge("FPN.LPFC_R", "FPN.PPC_R", 0.88, "Superior Longitudinal Fasciculus", PathwayModality::Association, false, None),
        edge("DAN.FEF_L", "DAN.IPS_L", 0.86, "Superior Longitudinal Fasciculus", PathwayModality::Association, false, None),
        edge("DAN.FEF_R", "DAN.IPS_R", 0.86, "Superior Longitudinal Fasciculus", PathwayModality::Association, false, None),
        edge("SAL.RPFC_L", "SAL.SMG_L", 0.74, "Superior Longitudinal Fasciculus", PathwayModality::Association, false, None),
        edge("SAL.RPFC_R", "SAL.SMG_R", 0.75, "Superior Longitudinal Fasciculus", PathwayModality::Association, false, None),
        edge("SMN.Lat_L", "DAN.IPS_L", 0.71, "Superior Longitudinal Fasciculus", PathwayModality::Association, false, None),
        edge("SMN.Lat_R", "DAN.IPS_R", 0.72, "Superior Longitudinal Fasciculus", PathwayModality::Association, false, None),
        // Uncinate Fasciculus (Frontotemporal Limbic)
        edge("SUBC.Amyg_L", "DMN.MPFC", 0.83, "Uncinate Fasciculus", PathwayModality::Association, false, None),
        edge("SUBC.Amyg_R", "DMN.MPFC", 0.81, "Uncinate Fasciculus", PathwayModality::Association, false, None),
        edge("SUBC.Amyg_L", "SAL.AIns_L", 0.78, "Uncinate Fasciculus", PathwayModality::Association, false, None),
        edge("SUBC.Amyg_R", "SAL.AIns_R", 0.77, "Uncinate Fasciculus", PathwayModality::Association, false, None),
        edge("SUBC.Hipp_L", "SUBC.Amyg_L", 0.85, "Uncinate Fasciculus", PathwayModality::Association, false, None),
        edge("SUBC.Amyg_L", "SAL.RPFC_L", 0.69, "Uncinate Fasciculus", PathwayModality::Association, false, None),
        // Inferior Fronto-Occipital Fasciculus (IFOF)
        edge("VIS.Occ", "VIS.Med", 0.93, "Inferior Fronto-Occipital Fasciculus", PathwayModality::Association, false, None),
        edge("VIS.Med", "VIS.Lat_L", 0.87, "Inferior Fronto-Occipital Fasciculus", PathwayModality::Association, false, None),
        edge("VIS.Med", "VIS.Lat_R", 0.87, "Inferior Fronto-Occipital Fasciculus", PathwayModality::Association, false, None),
        edge("VIS.Lat_L", "FPN.LPFC_L", 0.68, "Inferior Fronto-Occipital Fasciculus", PathwayModality::Association, false, None),
        edge("VIS.Lat_R", "FPN.LPFC_R", 0.67, "Inferior Fronto-Occipital Fasciculus", PathwayModality::Association, false, None),
        edge("VIS.Med", "DMN.PCC", 0.74, "Inferior Fronto-Occipital Fasciculus", PathwayModality::Association, false, None),
        // Corticospinal & Thalamocortical Projection Tracts
        edge("SUBC.Thal_L", "SMN.Lat_L", 0.86, "Corticospinal & Thalamocortical", PathwayModality::Projection, false, None),
        edge("SUBC.Thal_L", "SMN.Sup", 0.84, "Corticospinal & Thalamocortical", PathwayModality::Projection, false, None),
        edge("SUBC.Thal_L", "SAL.ACC", 0.76, "Corticospinal & Thalamocortical", PathwayModality::Projection, false, None),
        edge("SUBC.Thal_L", "FPN.LPFC_L", 0.75, "Corticospinal & Thalamocortical", PathwayModality::Projection, false, None),
        edge("SUBC.Caud_L", "SUBC.Pall_L", 0.88, "Corticospinal & Thalamocortical", PathwayModality::Projection, true, None),
        edge("SUBC.Pall_L", "SUBC.Thal_L", 0.85, "Corticospinal & Thalamocortical", PathwayModality::Projection, true, None),
        edge("FPN.LPFC_L", "SUBC.Caud_L", 0.79, "Corticospinal & Thalamocortical", PathwayModality::Projection, true, None),
        edge("SMN.Sup", "SMN.Lat_L", 0.83, "Corticospinal & Thalamocortical", PathwayModality::Projection, false, None),
        edge("SMN.Sup", "SMN.Lat_R", 0.83, "Corticospinal & Thalamocortical", PathwayModality::Projection, false, None),
        edge("SMN.Lat_L", "SUBC.Put_L", 0.86, "Corticospinal & Thalamocortical", PathwayModality::Projection, true, None),
        edge("SUBC.Put_L", "SUBC.Pall_L", 0.89, "Corticospinal & Thalamocortical", PathwayModality::Projection, true, None),
        edge("SMN.Lat_L", "SUBC.RN_L", 0.74, "Corticospinal & Thalamocortical", PathwayModality::Projection, true, None),
        edge("SUBC.Thal_L", "SUBC.RN_L", 0.72, "Corticospinal & Thalamocortical", PathwayModality::Projection, false, None),
        edge("FPN.LPFC_L", "BSTEM.Pons", 0.82, "Cerebellar Peduncles (CPCT / DRTT)", PathwayModality::Projection, true, None),
        edge("SMN.Lat_L", "BSTEM.Pons", 0.85, "Cerebellar Peduncles (CPCT / DRTT)", PathwayModality::Projection, true, None),
        edge("BSTEM.Pons", "CBL.Crus1_R", 0.88, "Cerebellar Peduncles (CPCT / DRTT)", PathwayModality::Projection, true, None),
        edge("BSTEM.Pons", "CBL.Crus1_L", 0.84, "Cerebellar Peduncles (CPCT / DRTT)", PathwayModality::Projection, true, None),
        edge("BSTEM.Pons", "CBL.LobVI_R", 0.87, "Cerebellar Peduncles (CPCT / DRTT)", PathwayModality::Projection, true, None),
        edge("CBL.Crus1_R", "SUBC.RN_L", 0.83, "Cerebellar Peduncles (CPCT / DRTT)", PathwayModality::Projection, true, None),
        edge("CBL.LobVI_R", "SUBC.Thal_L", 0.84, "Cerebellar Peduncles (CPCT / DRTT)", PathwayModality::Projection, true, None),
        edge("CBL.Crus1_L", "CBL.Crus1_R", 0.76, "Cerebellar Peduncles (CPCT / DRTT)", PathwayModality::Commissural, false, None),
        edge("SAL.ACC", "SAL.AIns_L", 0.84, "Cingulum Bundle", PathwayModality::Association, false, None),
        edge("SAL.ACC", "SAL.AIns_R", 0.83, "Cingulum Bundle", PathwayModality::Association, false, None),
    ];

    ConnectomeAtlas {
        id: "human_mni152".to_string(),
        title: "Human Macro-Connectome (MNI152 + Major Tracts)".to_string(),
        organism: "Homo sapiens".to_string(),
        coordinate_space: "MNI152 Stereotaxic Space (mm, +X Right, +Y Anterior, +Z Superior)".to_string(),
        scale: AtlasScale::MacroParcellation,
        reference_summary: "35 stereotaxic cortical and subcortical ROIs across 8 canonical intrinsic networks coupled by 7 major white-matter tract systems; contextualized with HCP-MMP1.0 (360 cortical + 14 subcortical parcels) and AAL3v1 (166 regions).".to_string(),
        primary_citation_url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC4990127/".to_string(),
        benchmark_neurons_or_parcels: "360 HCP-MMP1.0 cortical + 14 subcortical (166 AAL3v1)".to_string(),
        benchmark_synapses_or_tracts: "7 major tract systems (~200M callosal axons)".to_string(),
        regions,
        edges,
        ablated_regions: BTreeSet::new(),
        severed_tracts: BTreeSet::new(),
    }
}

/// 2. Mouse Mesoscale Connectome (Allen Mouse Brain CCFv3).
pub fn mouse_ccfv3_atlas() -> ConnectomeAtlas {
    let url_ccf = "https://alleninstitute.github.io/CCF-MAP/descriptions/mouse_ccf.html";
    let url_iso = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8026112/";
    let url_tc = "https://elifesciences.org/reviewed-preprints/93481v2/pdf";
    let url_bg = "https://pmc.ncbi.nlm.nih.gov/articles/PMC13156036/";

    let regions = vec![
        reg("MOB_L", "Main Olfactory Bulb (L)", "Olfactory", Hemisphere::Left, [-12.0, 46.0, -4.0], 4.6, "#fb923c", url_ccf, "CCFv3 main olfactory bulb annotated layers"),
        reg("PIR_L", "Piriform Cortex (L)", "Olfactory", Hemisphere::Left, [-28.0, 14.0, -18.0], 4.2, "#fb923c", url_iso, "Primary olfactory piriform cortex"),
        reg("FRP_L", "Frontal Pole (L)", "Prefrontal Isocortex", Hemisphere::Left, [-10.0, 34.0, 8.0], 3.8, "#facc15", url_tc, "Rostral frontal pole isocortex"),
        reg("ORB_L", "Orbital Area (L)", "Prefrontal Isocortex", Hemisphere::Left, [-14.0, 28.0, -4.0], 4.2, "#facc15", url_iso, "Orbitofrontal area (ORBm/ORBl)"),
        reg("ACA_L", "Anterior Cingulate Area (L)", "Medial Isocortex", Hemisphere::Left, [-5.0, 14.0, 16.0], 4.4, "#38bdf8", url_tc, "Dorsal/ventral anterior cingulate cortex"),
        reg("PL_L", "Prelimbic Area (L)", "Prefrontal Isocortex", Hemisphere::Left, [-4.0, 24.0, 6.0], 4.0, "#facc15", url_iso, "Medial prefrontal prelimbic cortex"),
        reg("ILA_L", "Infralimbic Area (L)", "Prefrontal Isocortex", Hemisphere::Left, [-4.0, 22.0, -4.0], 3.8, "#facc15", url_iso, "Ventromedial infralimbic cortex"),
        reg("MOs_L", "Secondary Motor Area (L)", "Sensorimotor Isocortex", Hemisphere::Left, [-12.0, 18.0, 20.0], 4.8, "#34d399", url_tc, "Premotor MOs (41.0% sensorimotor input from MOp)"),
        reg("MOs_R", "Secondary Motor Area (R)", "Sensorimotor Isocortex", Hemisphere::Right, [12.0, 18.0, 20.0], 4.8, "#34d399", url_tc, "Contralateral MOs homotypic callosal partner"),
        reg("MOp_L", "Primary Motor Area (L)", "Sensorimotor Isocortex", Hemisphere::Left, [-20.0, 6.0, 20.0], 5.2, "#34d399", url_tc, "Primary motor cortex MOp (63.2% sensory input from SSp)"),
        reg("MOp_R", "Primary Motor Area (R)", "Sensorimotor Isocortex", Hemisphere::Right, [20.0, 6.0, 20.0], 5.2, "#34d399", url_tc, "Contralateral primary motor cortex MOp"),
        reg("SSp_L", "Primary Somatosensory Area (L)", "Sensorimotor Isocortex", Hemisphere::Left, [-30.0, -4.0, 16.0], 5.4, "#34d399", url_tc, "Barrel/limb primary somatosensory cortex SSp"),
        reg("SSp_R", "Primary Somatosensory Area (R)", "Sensorimotor Isocortex", Hemisphere::Right, [30.0, -4.0, 16.0], 5.4, "#34d399", url_tc, "Contralateral primary somatosensory cortex SSp"),
        reg("SSs_L", "Supplemental Somatosensory (L)", "Sensorimotor Isocortex", Hemisphere::Left, [-36.0, -8.0, 6.0], 4.4, "#34d399", url_tc, "Lateral supplemental somatosensory area SSs"),
        reg("AI_L", "Agranular Insular Area (L)", "Medial Isocortex", Hemisphere::Left, [-32.0, 12.0, -6.0], 4.1, "#38bdf8", url_tc, "Insular viscerosensory/salience area"),
        reg("RSP_L", "Retrosplenial Area (L)", "Medial Isocortex", Hemisphere::Left, [-6.0, -24.0, 18.0], 4.6, "#38bdf8", url_iso, "Posterior medial retrosplenial navigation hub"),
        reg("AUD_L", "Auditory Areas (L)", "Sensory Isocortex", Hemisphere::Left, [-38.0, -22.0, 4.0], 4.4, "#a78bfa", url_tc, "Temporal auditory isocortical areas"),
        reg("VISp_L", "Primary Visual Area (L)", "Sensory Isocortex", Hemisphere::Left, [-26.0, -36.0, 14.0], 5.0, "#a78bfa", url_iso, "Posterior primary visual area VISp"),
        reg("VISl_L", "Lateral Visual Area (L)", "Sensory Isocortex", Hemisphere::Left, [-34.0, -34.0, 10.0], 4.2, "#a78bfa", url_iso, "Extrastriate lateral visual area VISl"),
        reg("CA1_L", "Ammon's Horn CA1 (L)", "Hippocampal Formation", Hemisphere::Left, [-20.0, -20.0, 6.0], 4.6, "#f472b6", url_ccf, "Dorsal hippocampal CA1 pyramidal layer"),
        reg("DG_L", "Dentate Gyrus (L)", "Hippocampal Formation", Hemisphere::Left, [-16.0, -22.0, 2.0], 4.2, "#f472b6", url_iso, "Hippocampal dentate gyrus granule layer"),
        reg("SUB_L", "Subiculum (L)", "Hippocampal Formation", Hemisphere::Left, [-18.0, -30.0, 4.0], 4.0, "#f472b6", url_iso, "Major hippocampal output subiculum"),
        reg("BLA_L", "Basolateral Amygdala (L)", "Cortical Subplate", Hemisphere::Left, [-28.0, -10.0, -16.0], 4.1, "#f87171", url_iso, "Basolateral amygdalar nucleus"),
        reg("CP_L", "Caudoputamen Striatum (L)", "Basal Ganglia", Hemisphere::Left, [-20.0, 6.0, 0.0], 5.2, "#2dd4bf", url_bg, "Dorsal striatum CP (5 CCFv3 multimodal subdivisions, D1/D2 MSNs)"),
        reg("GPe_L", "Globus Pallidus External (L)", "Basal Ganglia", Hemisphere::Left, [-22.0, -4.0, -6.0], 4.0, "#2dd4bf", url_bg, "External pallidum projecting within BG and sparsely to cortex"),
        reg("VAL_L", "Ventral Anterior-Lateral Thalamus (L)", "Thalamus", Hemisphere::Left, [-12.0, -10.0, -2.0], 4.3, "#fb7185", url_tc, "Motor thalamus VAL projecting strongly to MOp and MOs"),
        reg("VPM_L", "Ventral Posteromedial Thalamus (L)", "Thalamus", Hemisphere::Left, [-16.0, -14.0, -4.0], 4.3, "#fb7185", url_tc, "Somatosensory thalamus VPM projecting to SSp"),
        reg("PF_L", "Parafascicular Thalamic Nucleus (L)", "Thalamus", Hemisphere::Left, [-8.0, -18.0, -4.0], 3.9, "#fb7185", url_bg, "Intralaminar PF hub linking cortex to CP subdivisions"),
    ];

    let edges = vec![
        // Verified sensorimotor reciprocals & callosal projections (eLife 93481v2)
        edge("MOp_L", "MOs_L", 0.91, "Corticocortical Sensorimotor", PathwayModality::Association, true, None),
        edge("SSp_L", "MOp_L", 0.94, "Corticocortical Sensorimotor", PathwayModality::Association, true, None),
        edge("SSp_L", "MOs_L", 0.82, "Corticocortical Sensorimotor", PathwayModality::Association, true, None),
        edge("SSs_L", "SSp_L", 0.88, "Corticocortical Sensorimotor", PathwayModality::Association, true, None),
        edge("SSs_L", "MOp_L", 0.76, "Corticocortical Sensorimotor", PathwayModality::Association, true, None),
        edge("MOs_L", "MOp_L", 0.79, "Corticocortical Sensorimotor", PathwayModality::Association, true, None),
        edge("MOp_L", "MOp_R", 0.74, "Corpus Callosum (cc)", PathwayModality::Commissural, false, None),
        edge("MOs_L", "MOs_R", 0.75, "Corpus Callosum (cc)", PathwayModality::Commissural, false, None),
        edge("SSp_L", "SSp_R", 0.71, "Corpus Callosum (cc)", PathwayModality::Commissural, false, None),
        // Thalamocortical & Corticothalamic Loops
        edge("VAL_L", "MOp_L", 0.90, "Thalamocortical Radiation", PathwayModality::Projection, true, None),
        edge("VAL_L", "MOs_L", 0.88, "Thalamocortical Radiation", PathwayModality::Projection, true, None),
        edge("VPM_L", "SSp_L", 0.93, "Thalamocortical Radiation", PathwayModality::Projection, true, None),
        edge("PF_L", "CP_L", 0.89, "Thalamostriatal Pathway", PathwayModality::Projection, true, None),
        edge("ACA_L", "PF_L", 0.81, "Corticothalamic Projection", PathwayModality::Projection, true, None),
        // Basal Ganglia Loop (PMC13156036)
        edge("MOp_L", "CP_L", 0.87, "Corticostriatal Tract", PathwayModality::Projection, true, None),
        edge("SSp_L", "CP_L", 0.84, "Corticostriatal Tract", PathwayModality::Projection, true, None),
        edge("ACA_L", "CP_L", 0.78, "Corticostriatal Tract", PathwayModality::Projection, true, None),
        edge("ORB_L", "CP_L", 0.76, "Corticostriatal Tract", PathwayModality::Projection, true, None),
        edge("CP_L", "GPe_L", 0.88, "Striatopallidal Pathway", PathwayModality::Projection, true, None),
        edge("GPe_L", "VAL_L", 0.77, "Pallidothalamic Pathway", PathwayModality::Projection, true, None),
        edge("GPe_L", "MOp_L", 0.58, "Pallidocortical Sparse Projection", PathwayModality::Projection, true, None),
        // Medial, Visual, Hippocampal & Olfactory Pathways
        edge("VISp_L", "VISl_L", 0.90, "Visual Dorsal/Lateral Stream", PathwayModality::Association, true, None),
        edge("VISp_L", "RSP_L", 0.85, "Visual-Retrosplenial Stream", PathwayModality::Association, true, None),
        edge("RSP_L", "ACA_L", 0.86, "Cingulum / Medial Stream", PathwayModality::Association, false, None),
        edge("ACA_L", "PL_L", 0.83, "Medial Prefrontal Stream", PathwayModality::Association, false, None),
        edge("PL_L", "ILA_L", 0.84, "Medial Prefrontal Stream", PathwayModality::Association, false, None),
        edge("FRP_L", "ORB_L", 0.79, "Prefrontal Association", PathwayModality::Association, false, None),
        edge("ORB_L", "AI_L", 0.77, "Orbito-Insular Stream", PathwayModality::Association, false, None),
        edge("AUD_L", "VISl_L", 0.68, "Audiovisual Association", PathwayModality::Association, false, None),
        edge("DG_L", "CA1_L", 0.92, "Trisynaptic Hippocampal Loop", PathwayModality::Association, true, None),
        edge("CA1_L", "SUB_L", 0.91, "Trisynaptic Hippocampal Loop", PathwayModality::Association, true, None),
        edge("SUB_L", "RSP_L", 0.82, "Subiculo-Retrosplenial Tract", PathwayModality::Association, true, None),
        edge("CA1_L", "PL_L", 0.75, "Hippocampo-Prefrontal Pathway", PathwayModality::Projection, true, None),
        edge("BLA_L", "PL_L", 0.81, "Amygdalo-Prefrontal Pathway", PathwayModality::Association, false, None),
        edge("BLA_L", "ILA_L", 0.80, "Amygdalo-Prefrontal Pathway", PathwayModality::Association, false, None),
        edge("BLA_L", "CA1_L", 0.78, "Amygdalo-Hippocampal Pathway", PathwayModality::Association, false, None),
        edge("MOB_L", "PIR_L", 0.93, "Lateral Olfactory Tract (lot)", PathwayModality::Projection, true, None),
        edge("PIR_L", "ORB_L", 0.77, "Olfacto-Orbitofrontal Pathway", PathwayModality::Association, true, None),
        edge("PIR_L", "BLA_L", 0.79, "Olfacto-Amygdalar Pathway", PathwayModality::Association, true, None),
    ];

    ConnectomeAtlas {
        id: "mouse_ccfv3".to_string(),
        title: "Mouse Mesoscale Connectome (Allen Brain CCFv3)".to_string(),
        organism: "Mus musculus (C57BL/6J)".to_string(),
        coordinate_space: "Allen CCFv3 Centered Stereotaxic Space (scaled 10x mm, Anterior Commissure Origin)".to_string(),
        scale: AtlasScale::MesoscaleProjection,
        reference_summary: "28 CCFv3 structures spanning isocortex, hippocampal formation, olfactory bulb, basolateral amygdala, caudoputamen, pallidum, and thalamus (VAL, VPM, PF), grounded in the 10 µm 1,675-mouse reference space.".to_string(),
        primary_citation_url: url_ccf.to_string(),
        benchmark_neurons_or_parcels: "43 isocortical + 314 subcortical structures (10 µm voxels)".to_string(),
        benchmark_synapses_or_tracts: "81 annotated white-matter fiber tracts + AAV viral projections".to_string(),
        regions,
        edges,
        ablated_regions: BTreeSet::new(),
        severed_tracts: BTreeSet::new(),
    }
}

/// 3. Adult Drosophila Whole-Brain Synaptic Connectome (FlyWire FAFB v783).
pub fn drosophila_flywire_atlas() -> ConnectomeAtlas {
    let url_fw = "https://www.brainfacts.org/neuroscience-in-society/supporting-research/2024/researchers-create-first-adult-fruit-fly-brain-connectome-110724";
    let url_codex = "https://codex.flywire.ai/faq";
    let url_mb = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7909955/";
    let url_gal = "https://flywire.ai/gallery";

    let regions = vec![
        // Bilateral Optic Lobes
        reg("ME_L", "Medulla (L)", "Optic Lobe", Hemisphere::Left, [-58.0, -6.0, 4.0], 6.0, "#a78bfa", url_gal, "Left optic lobe medulla; highest cross-brain visual connectivity"),
        reg("ME_R", "Medulla (R)", "Optic Lobe", Hemisphere::Right, [58.0, -6.0, 4.0], 6.0, "#a78bfa", url_gal, "Right optic lobe medulla"),
        reg("LO_L", "Lobula Complex (L)", "Optic Lobe", Hemisphere::Left, [-44.0, -10.0, 2.0], 4.8, "#a78bfa", url_gal, "Left motion/feature-selective lobula neuropil"),
        reg("LO_R", "Lobula Complex (R)", "Optic Lobe", Hemisphere::Right, [44.0, -10.0, 2.0], 4.8, "#a78bfa", url_gal, "Right motion/feature-selective lobula neuropil"),
        reg("AOTU_L", "Anterior Optic Tubercle (L)", "Optic Lobe", Hemisphere::Left, [-28.0, 16.0, 14.0], 4.0, "#a78bfa", url_gal, "Sky-compass visual relay to central complex"),
        reg("AOTU_R", "Anterior Optic Tubercle (R)", "Optic Lobe", Hemisphere::Right, [28.0, 16.0, 14.0], 4.0, "#a78bfa", url_gal, "Right anterior optic tubercle"),
        // Olfactory & Associative Mushroom Body / Lateral Horn
        reg("AL_L", "Antennal Lobe (L)", "Olfactory System", Hemisphere::Left, [-16.0, 26.0, -10.0], 4.6, "#fb923c", url_mb, "51 olfactory glomeruli (~190 PN types)"),
        reg("AL_R", "Antennal Lobe (R)", "Olfactory System", Hemisphere::Right, [16.0, 26.0, -10.0], 4.6, "#fb923c", url_mb, "Right antennal lobe glomeruli"),
        reg("MB_CA_L", "Mushroom Body Calyx (L)", "Associative Memory", Hemisphere::Left, [-22.0, -18.0, 22.0], 4.5, "#facc15", url_mb, "Kenyon cell dendritic calyx receiving olfactory PN input"),
        reg("MB_CA_R", "Mushroom Body Calyx (R)", "Associative Memory", Hemisphere::Right, [22.0, -18.0, 22.0], 4.5, "#facc15", url_mb, "Right mushroom body calyx"),
        reg("MB_Lobes_L", "Mushroom Body Lobes (L)", "Associative Memory", Hemisphere::Left, [-18.0, 14.0, 12.0], 5.0, "#facc15", url_mb, "Kenyon axon compartments + DAN/MBON valence synapses"),
        reg("MB_Lobes_R", "Mushroom Body Lobes (R)", "Associative Memory", Hemisphere::Right, [18.0, 14.0, 12.0], 5.0, "#facc15", url_mb, "Right mushroom body lobes"),
        reg("LH_L", "Lateral Horn (L)", "Olfactory System", Hemisphere::Left, [-34.0, -4.0, 18.0], 4.4, "#fb923c", url_mb, "Innate olfactory valence & multimodal center"),
        reg("LH_R", "Lateral Horn (R)", "Olfactory System", Hemisphere::Right, [34.0, -4.0, 18.0], 4.4, "#fb923c", url_mb, "Right lateral horn"),
        // Central Complex (Navigation & Heading Ring Attractor)
        reg("PB", "Protocerebral Bridge", "Central Complex", Hemisphere::Midline, [0.0, -20.0, 20.0], 4.6, "#38bdf8", url_gal, "Posterior central complex heading bridge"),
        reg("FB", "Fan-Shaped Body", "Central Complex", Hemisphere::Midline, [0.0, -4.0, 12.0], 5.4, "#38bdf8", url_gal, "Layered vector navigation & goal neuropil"),
        reg("EB", "Ellipsoid Body", "Central Complex", Hemisphere::Midline, [0.0, 6.0, 4.0], 4.8, "#38bdf8", url_gal, "Toroidal head-direction ring-attractor neuropil"),
        reg("NO", "Noduli", "Central Complex", Hemisphere::Midline, [0.0, -6.0, -2.0], 3.8, "#38bdf8", url_gal, "Paired self-motion angular velocity neuropils"),
        // Mechanosensory, Premotor & Descending
        reg("AMMC_L", "Antennal Mechanosensory & Motor (L)", "Mechanosensory & Premotor", Hemisphere::Left, [-22.0, 10.0, -20.0], 4.3, "#34d399", url_gal, "Auditory courtship song & Johnston's organ hub"),
        reg("AVLP_L", "Anterior Ventrolateral Protocerebrum (L)", "Mechanosensory & Premotor", Hemisphere::Left, [-30.0, 8.0, -4.0], 4.5, "#34d399", url_gal, "Sensorimotor integration & auditory/visual relay"),
        reg("AVLP_R", "Anterior Ventrolateral Protocerebrum (R)", "Mechanosensory & Premotor", Hemisphere::Right, [30.0, 8.0, -4.0], 4.5, "#34d399", url_gal, "Right AVLP premotor hub"),
        reg("GNG", "Gnathal Ganglia", "Mechanosensory & Premotor", Hemisphere::Midline, [0.0, 12.0, -30.0], 5.6, "#34d399", url_gal, "Ventral gustatory & descending motor command hub"),
    ];

    let edges = vec![
        // Olfactory Antennal Lobe -> Mushroom Body & Lateral Horn (mALT)
        edge("AL_L", "MB_CA_L", 0.94, "Medial Antennoprotocerebral Tract (mALT)", PathwayModality::ChemicalSynapse, true, Some(420_000)),
        edge("AL_R", "MB_CA_R", 0.94, "Medial Antennoprotocerebral Tract (mALT)", PathwayModality::ChemicalSynapse, true, Some(418_000)),
        edge("AL_L", "LH_L", 0.89, "Medial Antennoprotocerebral Tract (mALT)", PathwayModality::ChemicalSynapse, true, Some(310_000)),
        edge("AL_R", "LH_R", 0.89, "Medial Antennoprotocerebral Tract (mALT)", PathwayModality::ChemicalSynapse, true, Some(308_000)),
        edge("MB_CA_L", "MB_Lobes_L", 0.96, "Kenyon Cell Pedunculus", PathwayModality::ChemicalSynapse, true, Some(1_250_000)),
        edge("MB_CA_R", "MB_Lobes_R", 0.96, "Kenyon Cell Pedunculus", PathwayModality::ChemicalSynapse, true, Some(1_240_000)),
        edge("MB_Lobes_L", "LH_L", 0.78, "MBON Convergent Valence Tract", PathwayModality::ChemicalSynapse, true, Some(185_000)),
        edge("MB_Lobes_R", "LH_R", 0.78, "MBON Convergent Valence Tract", PathwayModality::ChemicalSynapse, true, Some(182_000)),
        edge("MB_Lobes_L", "MB_Lobes_R", 0.74, "Interhemispheric MB Commissure", PathwayModality::Commissural, false, Some(140_000)),
        // Optic Lobe -> Central Brain & Contralateral Medulla
        edge("ME_L", "LO_L", 0.95, "Optic Chiasm Pathway", PathwayModality::ChemicalSynapse, true, Some(2_400_000)),
        edge("ME_R", "LO_R", 0.95, "Optic Chiasm Pathway", PathwayModality::ChemicalSynapse, true, Some(2_390_000)),
        edge("ME_L", "ME_R", 0.83, "Posterior Optic Commissure", PathwayModality::Commissural, false, Some(290_000)),
        edge("ME_L", "AOTU_L", 0.86, "Anterior Optic Tract", PathwayModality::ChemicalSynapse, true, Some(260_000)),
        edge("ME_R", "AOTU_R", 0.86, "Anterior Optic Tract", PathwayModality::ChemicalSynapse, true, Some(258_000)),
        edge("AOTU_L", "EB", 0.88, "Tuberculo-Bulbar Ring Pathway", PathwayModality::ChemicalSynapse, true, Some(215_000)),
        edge("AOTU_R", "EB", 0.88, "Tuberculo-Bulbar Ring Pathway", PathwayModality::ChemicalSynapse, true, Some(214_000)),
        edge("LO_L", "AVLP_L", 0.85, "Lobula Columnar Visuomotor Tract", PathwayModality::ChemicalSynapse, true, Some(340_000)),
        edge("LO_R", "AVLP_R", 0.85, "Lobula Columnar Visuomotor Tract", PathwayModality::ChemicalSynapse, true, Some(335_000)),
        // Central Complex Navigation Circuit
        edge("PB", "EB", 0.93, "Columnar E-PG / P-EN Compass Tract", PathwayModality::ChemicalSynapse, false, Some(480_000)),
        edge("PB", "FB", 0.91, "Protocerebral-Fan Columnar Tract", PathwayModality::ChemicalSynapse, true, Some(440_000)),
        edge("FB", "NO", 0.90, "Fan-Shaped to Noduli Tract", PathwayModality::ChemicalSynapse, true, Some(275_000)),
        edge("EB", "FB", 0.84, "Central Complex Ring-to-Fan Relay", PathwayModality::ChemicalSynapse, false, Some(230_000)),
        edge("MB_Lobes_L", "FB", 0.76, "MBON-to-Central-Complex Goal Pathway", PathwayModality::ChemicalSynapse, true, Some(150_000)),
        edge("MB_Lobes_R", "FB", 0.76, "MBON-to-Central-Complex Goal Pathway", PathwayModality::ChemicalSynapse, true, Some(148_000)),
        // Premotor & Descending Motor Output to Gnathal Ganglia
        edge("AMMC_L", "AVLP_L", 0.84, "Mechanosensory Wedge-AVLP Tract", PathwayModality::ChemicalSynapse, true, Some(195_000)),
        edge("AMMC_L", "GNG", 0.82, "Mechanosensory-Gnathal Tract", PathwayModality::ChemicalSynapse, true, Some(175_000)),
        edge("LH_L", "AVLP_L", 0.80, "Lateral Horn Premotor Output", PathwayModality::ChemicalSynapse, true, Some(165_000)),
        edge("LH_R", "AVLP_R", 0.80, "Lateral Horn Premotor Output", PathwayModality::ChemicalSynapse, true, Some(164_000)),
        edge("AVLP_L", "AVLP_R", 0.75, "Protocerebral Commissure", PathwayModality::Commissural, false, Some(130_000)),
        edge("AVLP_L", "GNG", 0.87, "Descending Premotor Tract", PathwayModality::Projection, true, Some(390_000)),
        edge("AVLP_R", "GNG", 0.87, "Descending Premotor Tract", PathwayModality::Projection, true, Some(388_000)),
        edge("FB", "GNG", 0.79, "Central Complex Steering Descending Tract", PathwayModality::Projection, true, Some(210_000)),
    ];

    ConnectomeAtlas {
        id: "drosophila_flywire".to_string(),
        title: "Drosophila Whole-Brain Connectome (FlyWire FAFB v783)".to_string(),
        organism: "Drosophila melanogaster (Adult Female)".to_string(),
        coordinate_space: "FAFB Central Brain & Optic Lobe Neuropil Space (scaled µm)".to_string(),
        scale: AtlasScale::MicroNeuropilSynaptic,
        reference_summary: "22 neuropil compartments across the adult fruit fly whole-brain EM connectome (139,255 proofread neurons, >50 million synaptic connections, >3,600 cell types across 78 neuropils).".to_string(),
        primary_citation_url: url_fw.to_string(),
        benchmark_neurons_or_parcels: "139,255 proofread neurons across 78 neuropils".to_string(),
        benchmark_synapses_or_tracts: ">50,000,000 EM-reconstructed chemical synapses (>3,600 cell types)".to_string(),
        regions: regions
            .into_iter()
            .map(|mut r| {
                if r.citation_url.is_empty() {
                    r.citation_url = url_codex.to_string();
                }
                r
            })
            .collect(),
        edges,
        ablated_regions: BTreeSet::new(),
        severed_tracts: BTreeSet::new(),
    }
}

/// 4. Complete Nematode Cellular Connectome (WormAtlas C. elegans Hermaphrodite).
pub fn celegans_wormatlas() -> ConnectomeAtlas {
    let url_wa = "https://www.wormatlas.org/neuronalwiring.html";
    let url_var = "https://arxiv.org/html/0907.2373v4";
    let url_cook = "https://pmc.ncbi.nlm.nih.gov/articles/PMC10983851/";
    let url_ov = "https://www.wormatlas.org/hermaphrodite/nervous/Neuroframeset.html";

    let regions = vec![
        // Sensory Neurons (Anterior Amphid, Gentle/Harsh Touch, Posterior Touch)
        reg("ASHL", "ASHL Amphid Nociceptor", "Sensory Neurons", Hemisphere::Left, [-8.0, 52.0, 6.0], 3.8, "#fb923c", url_ov, "Left amphid polymodal nociceptive/avoidance sensory neuron"),
        reg("ASHR", "ASHR Amphid Nociceptor", "Sensory Neurons", Hemisphere::Right, [8.0, 52.0, 6.0], 3.8, "#fb923c", url_ov, "Right amphid polymodal nociceptive sensory neuron"),
        reg("ALML", "ALML Anterior Touch", "Sensory Neurons", Hemisphere::Left, [-14.0, 18.0, 10.0], 3.8, "#fb923c", url_var, "Left anterior lateral microtubule gentle-touch mechanoreceptor"),
        reg("ALMR", "ALMR Anterior Touch", "Sensory Neurons", Hemisphere::Right, [14.0, 18.0, 10.0], 3.8, "#fb923c", url_var, "Right anterior lateral microtubule gentle-touch mechanoreceptor"),
        reg("AVM", "AVM Ventral Touch", "Sensory Neurons", Hemisphere::Midline, [0.0, 22.0, -8.0], 3.8, "#fb923c", url_var, "Anterior ventral microtubule gentle-touch receptor"),
        reg("PLML", "PLML Posterior Touch", "Sensory Neurons", Hemisphere::Left, [-12.0, -56.0, -4.0], 3.8, "#fb923c", url_var, "Left posterior lateral gentle-touch mechanoreceptor"),
        reg("PLMR", "PLMR Posterior Touch", "Sensory Neurons", Hemisphere::Right, [12.0, -56.0, -4.0], 3.8, "#fb923c", url_var, "Right posterior lateral gentle-touch mechanoreceptor"),
        reg("PVDL", "PVDL Harsh Touch Polymodal", "Sensory Neurons", Hemisphere::Left, [-16.0, -24.0, 6.0], 3.6, "#fb923c", url_cook, "Left multidendritic harsh-touch & cold nociceptor"),
        // Rich-Club Command Interneurons (AVA, AVB, AVD, AVE, PVC; mean combined indegree 48.5)
        reg("AVAL", "AVAL Backward Command Hub", "Rich-Club Command", Hemisphere::Left, [-6.0, 36.0, 4.0], 5.2, "#38bdf8", url_var, "Rich-club command interneuron coordinating backward locomotion"),
        reg("AVAR", "AVAR Backward Command Hub", "Rich-Club Command", Hemisphere::Right, [6.0, 36.0, 4.0], 5.2, "#38bdf8", url_var, "Right rich-club backward locomotion command hub"),
        reg("AVBL", "AVBL Forward Command Hub", "Rich-Club Command", Hemisphere::Left, [-7.0, 32.0, -2.0], 5.1, "#38bdf8", url_var, "Rich-club command interneuron coordinating forward locomotion"),
        reg("AVBR", "AVBR Forward Command Hub", "Rich-Club Command", Hemisphere::Right, [7.0, 32.0, -2.0], 5.1, "#38bdf8", url_var, "Right rich-club forward locomotion command hub"),
        reg("AVDL", "AVDL Reversal Command", "Rich-Club Command", Hemisphere::Left, [-9.0, 28.0, 8.0], 4.5, "#38bdf8", url_cook, "Anterior touch-evoked backward command interneuron"),
        reg("AVDR", "AVDR Reversal Command", "Rich-Club Command", Hemisphere::Right, [9.0, 28.0, 8.0], 4.5, "#38bdf8", url_cook, "Right reversal command interneuron"),
        reg("AVEL", "AVEL Head/Reversal Command", "Rich-Club Command", Hemisphere::Left, [-5.0, 40.0, 0.0], 4.4, "#38bdf8", url_ov, "Ventral cord anterior reversal command interneuron"),
        reg("PVCL", "PVCL Forward Command (Tail)", "Rich-Club Command", Hemisphere::Left, [-8.0, -48.0, -2.0], 4.7, "#38bdf8", url_cook, "Lumbar ganglion rich-club forward locomotion command interneuron"),
        reg("PVCR", "PVCR Forward Command (Tail)", "Rich-Club Command", Hemisphere::Right, [8.0, -48.0, -2.0], 4.7, "#38bdf8", url_cook, "Right lumbar forward command interneuron"),
        // Nerve-Ring & Somatic-Pharyngeal Bridge Interneurons
        reg("RIBL", "RIBL Nerve-Ring Integrator", "Nerve-Ring & Bridge", Hemisphere::Left, [-10.0, 42.0, -4.0], 4.3, "#facc15", url_var, "High-centrality nerve-ring amphid integration interneuron"),
        reg("RIPL", "RIPL Pharyngeal Bridge", "Nerve-Ring & Bridge", Hemisphere::Left, [-4.0, 48.0, 2.0], 4.0, "#facc15", url_ov, "Sole somatic-to-pharyngeal nervous system bridge interneuron (L)"),
        reg("RIPR", "RIPR Pharyngeal Bridge", "Nerve-Ring & Bridge", Hemisphere::Right, [4.0, 48.0, 2.0], 4.0, "#facc15", url_ov, "Sole somatic-to-pharyngeal nervous system bridge interneuron (R)"),
        // Ventral Cord, Enteric & Pharyngeal Motor Neurons
        reg("DA01", "DA01 Dorsal Backward Motor", "Motor & Pharyngeal", Hemisphere::Midline, [0.0, 8.0, 12.0], 4.1, "#34d399", url_cook, "A-class cholinergic backward dorsal motor neuron"),
        reg("VA01", "VA01 Ventral Backward Motor", "Motor & Pharyngeal", Hemisphere::Midline, [0.0, 6.0, -12.0], 4.1, "#34d399", url_cook, "A-class cholinergic backward ventral motor neuron"),
        reg("DB01", "DB01 Dorsal Forward Motor", "Motor & Pharyngeal", Hemisphere::Midline, [0.0, -6.0, 12.0], 4.1, "#34d399", url_cook, "B-class cholinergic forward dorsal motor neuron"),
        reg("VB01", "VB01 Ventral Forward Motor", "Motor & Pharyngeal", Hemisphere::Midline, [0.0, -10.0, -12.0], 4.1, "#34d399", url_cook, "B-class cholinergic forward ventral motor neuron"),
        reg("DD02", "DD02 GABA Cross-Inhibitor", "Motor & Pharyngeal", Hemisphere::Midline, [0.0, -18.0, 0.0], 3.9, "#a78bfa", url_cook, "D-class inhibitory GABAergic dorsoventral cross-inhibitor"),
        reg("M4_PHAR", "M4 Pharyngeal Peristalsis Motor", "Motor & Pharyngeal", Hemisphere::Midline, [0.0, 60.0, 0.0], 4.2, "#f472b6", url_ov, "Pharyngeal nervous system (20-neuron subsystem) isthmus motor hub"),
    ];

    let edges = vec![
        // Sensory -> Command Chemical Synapses & Gap Junctions
        edge("ASHL", "AVAL", 0.89, "Amphid Avoidance Chemical", PathwayModality::ChemicalSynapse, true, Some(18)),
        edge("ASHR", "AVAR", 0.89, "Amphid Avoidance Chemical", PathwayModality::ChemicalSynapse, true, Some(17)),
        edge("ASHL", "AVDL", 0.82, "Amphid Avoidance Chemical", PathwayModality::ChemicalSynapse, true, Some(12)),
        edge("ASHR", "AVDR", 0.82, "Amphid Avoidance Chemical", PathwayModality::ChemicalSynapse, true, Some(12)),
        edge("ASHL", "RIBL", 0.78, "Amphid Nerve-Ring Chemical", PathwayModality::ChemicalSynapse, true, Some(10)),
        edge("ALML", "AVDL", 0.86, "Anterior Touch Reversal Chemical", PathwayModality::ChemicalSynapse, true, Some(15)),
        edge("ALMR", "AVDR", 0.86, "Anterior Touch Reversal Chemical", PathwayModality::ChemicalSynapse, true, Some(15)),
        edge("AVM", "AVBL", 0.74, "Touch-Command Gap Junction", PathwayModality::ElectricalGapJunction, false, Some(8)),
        edge("ALML", "PVCL", 0.79, "Anterior Touch Chemical", PathwayModality::ChemicalSynapse, true, Some(11)),
        edge("PLML", "PVCL", 0.93, "Posterior Touch Gap Junction", PathwayModality::ElectricalGapJunction, false, Some(22)),
        edge("PLMR", "PVCR", 0.93, "Posterior Touch Gap Junction", PathwayModality::ElectricalGapJunction, false, Some(22)),
        edge("PLML", "AVAL", 0.72, "Posterior Touch Chemical", PathwayModality::ChemicalSynapse, true, Some(7)),
        edge("PVDL", "AVAL", 0.84, "Harsh Touch Nociceptive Chemical", PathwayModality::ChemicalSynapse, true, Some(14)),
        edge("PVDL", "PVCL", 0.83, "Harsh Touch Nociceptive Chemical", PathwayModality::ChemicalSynapse, true, Some(13)),
        // Rich-Club Command Interneuron Core
        edge("AVAL", "AVAR", 0.95, "Rich-Club Command Gap Junction", PathwayModality::ElectricalGapJunction, false, Some(31)),
        edge("AVBL", "AVBR", 0.94, "Rich-Club Command Gap Junction", PathwayModality::ElectricalGapJunction, false, Some(29)),
        edge("AVDL", "AVAL", 0.91, "Rich-Club Command Chemical", PathwayModality::ChemicalSynapse, true, Some(26)),
        edge("AVDR", "AVAR", 0.91, "Rich-Club Command Chemical", PathwayModality::ChemicalSynapse, true, Some(25)),
        edge("AVEL", "AVAL", 0.87, "Rich-Club Command Chemical", PathwayModality::ChemicalSynapse, true, Some(19)),
        edge("PVCL", "AVBL", 0.90, "Rich-Club Command Chemical", PathwayModality::ChemicalSynapse, true, Some(24)),
        edge("PVCR", "AVBR", 0.90, "Rich-Club Command Chemical", PathwayModality::ChemicalSynapse, true, Some(24)),
        edge("PVCL", "AVAL", 0.81, "Rich-Club Cross-Command Chemical", PathwayModality::ChemicalSynapse, true, Some(14)),
        edge("RIBL", "AVBL", 0.85, "Nerve-Ring to Forward Command", PathwayModality::ElectricalGapJunction, false, Some(16)),
        edge("RIBL", "AVEL", 0.83, "Nerve-Ring to Reversal Command", PathwayModality::ChemicalSynapse, true, Some(15)),
        // Somatic-Pharyngeal Bridge via RIPL/RIPR
        edge("RIBL", "RIPL", 0.77, "Somatic-Pharyngeal Bridge", PathwayModality::ChemicalSynapse, true, Some(9)),
        edge("RIPL", "RIPR", 0.80, "Somatic-Pharyngeal Bridge", PathwayModality::ElectricalGapJunction, false, Some(10)),
        edge("RIPL", "M4_PHAR", 0.85, "Somatic-Pharyngeal Bridge", PathwayModality::ChemicalSynapse, true, Some(12)),
        edge("RIPR", "M4_PHAR", 0.85, "Somatic-Pharyngeal Bridge", PathwayModality::ChemicalSynapse, true, Some(12)),
        // Ventral Cord Motor Execution (62% chemical & 55% electrical input from command interneurons)
        edge("AVAL", "DA01", 0.94, "Ventral Cord Backward Drive", PathwayModality::ElectricalGapJunction, false, Some(28)),
        edge("AVAR", "VA01", 0.93, "Ventral Cord Backward Drive", PathwayModality::ChemicalSynapse, true, Some(27)),
        edge("AVEL", "DA01", 0.86, "Ventral Cord Backward Drive", PathwayModality::ChemicalSynapse, true, Some(16)),
        edge("AVBL", "DB01", 0.94, "Ventral Cord Forward Drive", PathwayModality::ElectricalGapJunction, false, Some(29)),
        edge("AVBR", "VB01", 0.93, "Ventral Cord Forward Drive", PathwayModality::ElectricalGapJunction, false, Some(28)),
        edge("DA01", "DD02", 0.84, "Dorsoventral GABA Cross-Inhibition", PathwayModality::ChemicalSynapse, true, Some(15)),
        edge("VB01", "DD02", 0.85, "Dorsoventral GABA Cross-Inhibition", PathwayModality::ChemicalSynapse, true, Some(16)),
    ];

    ConnectomeAtlas {
        id: "celegans_wormatlas".to_string(),
        title: "C. elegans Cellular Connectome (WormAtlas Hermaphrodite)".to_string(),
        organism: "Caenorhabditis elegans (Adult Hermaphrodite)".to_string(),
        coordinate_space: "Body-Axis Ganglion/Ventral-Cord Space (Anterior +Y Head/Pharynx, Posterior -Y Tail)".to_string(),
        scale: AtlasScale::CellularCompleteConnectome,
        reference_summary: "26 identifiable sensory, rich-club command (AVA, AVB, AVD, AVE, PVC), RIP pharyngeal bridge, and ventral cord motor neurons from the complete 302-neuron hermaphrodite wiring diagram.".to_string(),
        primary_citation_url: url_wa.to_string(),
        benchmark_neurons_or_parcels: "302 total neurons (282 somatic + 20 pharyngeal via RIPL/R)".to_string(),
        benchmark_synapses_or_tracts: "6,393 chemical synapses + 890 gap junctions + 1,410 NMJs".to_string(),
        regions,
        edges,
        ablated_regions: BTreeSet::new(),
        severed_tracts: BTreeSet::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_four_presets_compute_valid_topology_and_paths() {
        let presets = all_presets();
        assert_eq!(presets.len(), 4);
        for mut atlas in presets {
            let summary = atlas.compute_summary_metrics(0.0);
            assert!(summary.active_nodes >= 20);
            assert!(summary.active_edges >= 30);
            assert!(summary.global_efficiency > 0.0);
            assert_eq!(summary.connected_components, 1);

            let first = atlas.regions[0].id.clone();
            let last = atlas.regions[atlas.regions.len() - 1].id.clone();
            let path = atlas.shortest_signaling_path(&first, &last, 0.0);
            let p = path.expect("Expected intact signaling path");
            assert!(p.conduction_latency_ms > 0.0);
            assert!(summary.mean_communicability > 0.0);

            let first_tract = atlas.unique_tracts()[0].clone();
            atlas.toggle_tract_sever(&first_tract);
            let delta = atlas.evaluate_lesion_impact(0.0);
            assert!(delta.lesioned_efficiency <= delta.intact_efficiency + 1e-9);
        }
    }
}
