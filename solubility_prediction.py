#Import the required packages for this code.
from rdkit import Chem
from rdkit.Chem import Draw
import collections
import textwrap
import pandas as pd
import numpy as np
import base64
import re
from io import BytesIO
import plotly.graph_objects as go
import plotly.io as pio
import webbrowser
from itertools import combinations
import json

# --- 1. Corrected Group Dictionary ---
first_order_groups = {
    # --- Complexity 5 ---
    "CON(CH_{3})_{2}": "C(=O)N(C)C",

    # --- Complexity 4 ---
    "CH_{3} CO-": "[CH3][C](=O)[#6]",
    "CH_{2} CO": "[CH2][C](=O)[#6]",
    "CH_{3} COO": "[CH3][C](=O)[O][!#6]",
    "CH_{2} COO": "[CH2][C](=O)[O][!#1]",
    "C_{2} H_{5} O_{2}": "[CH2][CH2][O][O]",
    "CCl_{3}": "C(Cl)(Cl)(Cl)",
    "CF_{3}": "C(F)(F)(F)",
    "CH_{2} NO_{2}": "[CH2][N+](=O)[O-]",
    "CHNO": "[CH][N+](=O)[O-]",
    "ACNO_{2}": "[c][N+](=O)[O-]",

    # --- Complexity 3 ---
    "CHO(aldehydes)": "[CH](=O)[#6]",
    "COOH": "C(=O)[OH]",
    "HCOO": "[C;H1](=O)[O]",
    "COO": "[C;!H1](=O)[O]",
    "CONH_{2}": "C(=O)[NH2]",
    "CHCl_{2}": "[CH](Cl)(Cl)",
    "CCl_{2}": "[C;H0](Cl)(Cl)",
    "Cl-(C=C)": "[C](=C)[Cl]",
    "CH_{2} CN": "[CH2]C#N",
    "CF_{2}": "[C;H1](F)(F)",
    "CH_{2}=C=C<": "[CH2]=C=C",
    "O=C=N-": "N=C=O",
    "SO_{2}": "S(=O)(=O)",

    # --- Complexity 2 ---
    ">C=N-": ["[C;H0]=[N]", "[c;H0]:[n]"],
    "-CH=CH-": "[C;H1;!a]=[C;H1;!a]", # Corrected: Specifically targets non-aromatic alkene groups
    "CH_{3} N": "[CH3]-[N;X3]",
    "=CH-CH_{2}": "[CH2]=[CH]",
    "CH≡C-": "[CH]#C",
    "C≡C": "C#C",
    "ACCH_{3}": "[c][CH3]",
    "ACCH_{2}-": "[c][CH2]",
    "OH (alcohol)": "[OH1;X2;!$(O-c);!$(O-C=O)]", # UPDATED: Non-greedy alcohol O
    "OH (phenol)": "[OH1;X2;$(O-c)]",           # UPDATED: Non-greedy phenol O
    "CH_{2} O(cyclic)": "[O;R][CH2;R]",
    "CH_{3} O": "[CH3][O;X2]",
    "CH_{2} O": "[CH2][O;X2]",
    "CHO(ethers)": "[CH1;X4][O;X2]",
    "CH_{2} NH_{2}": "[CH2][NH2]",
    "CHNH_{2}": "[CH][NH2]",
    "CH_{3} NH": "[CH3][NH]",
    "CH_{2} NH": "[CH2][NH]",
    "CH_{2} N": "[CH2][N]",
    "ACNH_{2}": "[c][NH2]",
    "CH_{2} SH": "[CH2][SH]",
    "CH_{3} S": "[CH3][S]",
    "CH_{2} S": "[CH2][S]",
    "CH_{2} Cl": "[CH2][Cl]",
    "CHCl": "[CH][Cl]",
    "CCl": "[C;H0][Cl]",
    "ACCl": "[c][Cl]",
    "ACF": "[c][F]",
    "CF": "[CH][F]",
    "CN": "C#N",
    ">C=S": "C=S",
    ">C=0": "C=O",

    # --- Complexity 1 (most general) ---
    "NH": ["[NH1;!n]", "[nH]"],
    "NH2": ["[NH2;!n]", "[nH2]"],
    "pyridine_N": "[n;X2;H0]",
    ">C<": "[C;H0]",
    "-CH_{3}": "[CH3;X4;H3]",
    "-CH_{2}": "[CH2;X4;H2]",
    "- CH<": "[CH1;X4;H1]",
    "ACH": "[cH]",
    "AC": "[c;H0]",
    "I": "[I]",
    "Br": "[Br]",
    "F": "[F]",
    "SH": "[SH]",
    "S": "[S]",
    "N": "[N;H0;!n;!+]",
}

# --- 2. Pattern Preparation Code ---
group_patterns = {}
for name, smarts_val in first_order_groups.items():
    if isinstance(smarts_val, str):
        patterns = [Chem.MolFromSmarts(smarts_val)]
    else:
        patterns = [Chem.MolFromSmarts(s) for s in smarts_val]
    group_patterns[name] = [p for p in patterns if p is not None]

def get_pattern_complexity(patterns):
    if not patterns:
        return (0, 0)
    p = max(patterns, key=lambda m: m.GetNumAtoms() if m else 0)
    return (p.GetNumAtoms(), p.GetNumBonds()) if p else (0, 0)

sorted_group_patterns = sorted(
    group_patterns.items(),
    key=lambda item: get_pattern_complexity(item[1]),
    reverse=True
)

def count_first_order_groups(smiles):
    """
    Decomposes a molecule into functional groups. If any atoms are un-matched,
    it provides details and prompts the user to classify them.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        print(f"Error: Could not parse SMILES string: {smiles}")
        return {}

    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION)
    except Exception as e:
        print(f"Error during manual sanitization for SMILES '{smiles}': {e}")
        return {}

    mol = Chem.AddHs(mol)
    matched_atoms = set()
    group_counts = collections.defaultdict(int)

    for name, patterns in sorted_group_patterns:
        all_matches_for_group = []
        for pattern in patterns:
            all_matches_for_group.extend(mol.GetSubstructMatches(pattern, uniquify=True))

        for match in all_matches_for_group:
            # This check now correctly ignores hydrogens and dummy atoms.
            heavy_atoms = [idx for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1]
            if not any(atom_idx in matched_atoms for atom_idx in heavy_atoms):
                matched_atoms.update(heavy_atoms)
                group_counts[name] += 1

    # --- Validation and Interactive Prompt Step ---
    # FIXED: This now ignores hydrogens (atomic #1) and polymer dummy atoms [*] (atomic #0).
    all_heavy_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1}
    unmatched_indices = all_heavy_atoms.difference(matched_atoms)

    if unmatched_indices:
        print("\n" + "="*50)
        print("🚨 VALIDATION FAILED: Some atoms could not be matched automatically.")
        print("="*50)

        available_groups = list(first_order_groups.keys())

        for idx in sorted(list(unmatched_indices)):
            atom = mol.GetAtomWithIdx(idx)
            neighbors = [n.GetSymbol() + str(n.GetIdx()) for n in atom.GetNeighbors()]
            details = f"Atom {atom.GetIdx()} ({atom.GetSymbol()}) connected to {neighbors}"

            print(f"\n--- Unmatched Atom Found ---")
            print(f"Details: {details}")
            print("Please classify this atom by choosing from the list of available groups:")

            print(textwrap.fill(", ".join(available_groups), width=80))

            while True:
                chosen_group = input(f"Enter functional group name for Atom {idx}: ").strip()
                if chosen_group in available_groups:
                    group_counts[chosen_group] += 1
                    print(f"-> Classified Atom {idx} as '{chosen_group}'.")
                    break
                else:
                    print(f"❌ Error: '{chosen_group}' is not a valid key. Please choose from the list above.")

        print("\n" + "="*50)
        print("✅ All atoms have now been classified.")
        print("="*50 + "\n")

    return dict(group_counts)

def show_structure(smiles, save=True):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES string.")
        return None

    img = Draw.MolToImage(mol)
    buf = BytesIO()
    img.save(buf, format="PNG")
    img_base64 = base64.b64encode(buf.getvalue()).decode("utf-8")
    return img_base64


# Read the CSV file
df = pd.read_csv('Hansen_solubility_first_order_contributions.csv')

# Optional: clean whitespace and standardize column names
df.columns = df.columns.str.strip().str.lower()
df['first-order groups'] = df['first-order groups'].str.strip()

# Set index to group name for easy lookup
group_contributions = df.set_index('first-order groups')[['delta_d', 'delta_p', 'delta_h']]

# Load the Hansen contributions table once, treating N/A as 0
hansen_df = pd.read_csv(
    'Hansen_solubility_first_order_contributions.csv',
    na_values=["N/A", "nan"]
).fillna(0)

# Convert the "First-order groups" column to match the functional group keys you use
hansen_df.set_index("First-order groups", inplace=True)

def calculate_hansen_from_smiles(smiles):
    group_counts = count_first_order_groups(smiles)
    print("group_counts:", group_counts)  # Debug output

    delta_d = delta_p = delta_h = 0.0
    group_contributions = {}

    for group, count in group_counts.items():
        if count > 0:
            if group in hansen_df.index:
                val_d = hansen_df.loc[group, "delta_d"]
                val_p = hansen_df.loc[group, "delta_p"]
                val_h = hansen_df.loc[group, "delta_h"]
                print(f"Adding {count} * {group} values: delta_d={val_d}, delta_p={val_p}, delta_h={val_h}")
                contribution = {
                    "count": count,
                    "delta_d": val_d * count,
                    "delta_p": val_p * count,
                    "delta_h": val_h * count
                }
                group_contributions[group] = contribution
                delta_d += contribution["delta_d"]
                delta_p += contribution["delta_p"]
                delta_h += contribution["delta_h"]
            else:
                print(f"Group '{group}' not found in hansen_df index")

    print("Totals:", delta_d, delta_p, delta_h)
    return delta_d, delta_p, delta_h, group_contributions

def calculate_hansen_for_copolymer(smiles_percentage_dictionary):
    """
    Calculates the Hansen Solubility Parameters for a copolymer based on the
    mole fractions of its constituent monomers.
    """
    total_delta_d = 0.0
    total_delta_p = 0.0
    total_delta_h = 0.0
    monomer_details = {}

    for smiles, percentage in smiles_percentage_dictionary.items():
        mole_fraction = percentage / 100.0

        # Calculate the contributions for each individual monomer
        sigma_d, sigma_p, sigma_h, group_contribs = calculate_hansen_from_smiles(smiles)

        # Store the details for the final report
        monomer_details[smiles] = {
            'mole_fraction': mole_fraction,
            'group_contributions': group_contribs
        }

        # Add the weighted contribution of this monomer to the total
        total_delta_d += sigma_d * mole_fraction
        total_delta_p += sigma_p * mole_fraction
        total_delta_h += sigma_h * mole_fraction

    # Apply the final calibration constants to the summed, weighted contributions
    mw = 41.053 # This is a constant from the original calculation method
    final_delta_d = total_delta_d + (0.0 * mw) + 17.3231
    final_delta_p = total_delta_p + (0 * mw) + 7.3548
    final_delta_h = total_delta_h + (0.0 * mw) + 7.9793

    return final_delta_d, final_delta_p, final_delta_h, monomer_details


def plot_and_calculate_values(smiles):

    """
    Performs a full analysis for a single molecule, replicating the advanced
    copolymer report layout and functionality with the specific required calculations.
    """
    # --- This list can be customized as needed ---
    common_solvent_list = [
        'acetic acid', 'acetone', 'acetonitrile', 'benzene', '1-butanol', '2-butanol',
        '2-butanone', 't-butyl alcohol', 'carbon tetrachloride', 'chlorobenzene', 'chloroform',
        'cyclohexane', '1,2-dichloroethane', 'diethylene glycol', 'diethyl ether', 'diglyme',
        '1,2-dimethoxyethane', 'dimethylformamide', 'dimethyl sulfoxide', '1,4-dioxane',
        'ethanol', 'ethyl acetate', 'ethylene glycol', 'glycerin', 'heptane',
        'hexamethylphosphoramide', 'hexane', 'methanol', 'methyl t-butyl ether',
        'methylene chloride', 'n-methyl-2-pyrrolidinone', 'nitromethane', 'pentane',
        '1-propanol', '2-propanol', 'pyridine', 'tetrahydrofuran', 'toluene',
        'triethyl amine', 'water', 'o-xylene', 'm-xylene', 'p-xylene'
    ]

    # --- Part 1: Hansen Parameter Calculation for a Single Molecule ---
    try:
        final_group_counts = count_first_order_groups(smiles)
        if not final_group_counts:
              print("Could not determine functional groups. Aborting.")
              return
    except Exception as e:
        print(f"Error during group counting for SMILES '{smiles}': {e}")
        return

    delta_d = delta_p = delta_h = 0.0
    group_contributions = {}

    for group, count in final_group_counts.items():
        if count > 0 and group in hansen_df.index:
            contribution = {
                "count": count,
                "delta_d": hansen_df.loc[group, "delta_d"] * count,
                "delta_p": hansen_df.loc[group, "delta_p"] * count,
                "delta_h": hansen_df.loc[group, "delta_h"] * count
            }
            group_contributions[group] = contribution
            delta_d += contribution["delta_d"]
            delta_p += contribution["delta_p"]
            delta_h += contribution["delta_h"]

    # --- CALCULATION RESTORED ---
    # The required hardcoded values are now back in the calculation.
    mw = 41.053
    final_delta_d = delta_d + (0.0 * mw) + 17.3231
    final_delta_p = delta_p + (0 * mw) + 7.3548
    final_delta_h = delta_h + (0.0 * mw) + 7.9793

    print(f"Final Hansen Parameters (δD, δP, δH): ({final_delta_d:.2f}, {final_delta_p:.2f}, {final_delta_h:.2f})")

    # --- Part 2: Data Preparation for Report ---
    try:
        solvent_df = pd.read_csv('Solvent_values_hansen_solubility.csv')
    except FileNotFoundError:
        print("Error: 'Solvent_values_hansen_solubility.csv' not found. Aborting.")
        return

    solvent_df['is_common'] = solvent_df['Solvent'].str.lower().isin([s.lower() for s in common_solvent_list])

    def hansen_distance(row):
        return np.sqrt(4 * (final_delta_d - row["δD"])**2 + (final_delta_p - row["δP"])**2 + (final_delta_h - row["δH"])**2)
    solvent_df["Ra"] = solvent_df.apply(hansen_distance, axis=1)

    # Filter for single solvent matching tables
    all_matches = solvent_df.sort_values("Ra")[["Solvent", "Ra", "δD", "δP", "δH", "is_common"]]
    common_solvents_only_df = solvent_df[solvent_df['is_common']].copy().reset_index(drop=True)
    best_matches = common_solvents_only_df.nsmallest(5, "Ra")[["Solvent", "Ra", "δD", "δP", "δH"]]
    worst_matches = common_solvents_only_df.nlargest(5, "Ra")[["Solvent", "Ra", "δD", "δP", "δH"]]

    # Helper function for two-component optimization
    def get_optimized_mixtures_html(solvents_to_use_df, title, is_detailed_view=False):
        print(f"Performing optimization for: {title}...")
        mixture_results = []
        if len(solvents_to_use_df) < 2: return f"<h3>{title}</h3><p>Not enough solvents to create mixtures.</p>"
        for (idx1, idx2) in combinations(solvents_to_use_df.index, 2):
            s1_data, s2_data = solvents_to_use_df.loc[idx1], solvents_to_use_df.loc[idx2]
            for vol_frac1 in np.arange(0.05, 1.0, 0.05):
                mix_dD = (s1_data['δD'] * vol_frac1) + (s2_data['δD'] * (1.0-vol_frac1)); mix_dP = (s1_data['δP'] * vol_frac1) + (s2_data['δP'] * (1.0-vol_frac1)); mix_dH = (s1_data['δH'] * vol_frac1) + (s2_data['δH'] * (1.0-vol_frac1))
                ra_mix = np.sqrt(4*(final_delta_d-mix_dD)**2 + (final_delta_p-mix_dP)**2 + (final_delta_h-mix_dH)**2)
                mixture_results.append({"Solvent 1": s1_data['Solvent'], "Solvent 2": s2_data['Solvent'], "Vol Frac 1": vol_frac1, "Ra": ra_mix})
        if not mixture_results: return f"<h3>{title}</h3><p>No mixtures could be generated.</p>"
        mixtures_df = pd.DataFrame(mixture_results)
        best_pairs_df = mixtures_df.loc[mixtures_df.groupby(['Solvent 1', 'Solvent 2'])['Ra'].idxmin()]
        top_20_best_pairs = best_pairs_df.nsmallest(20, 'Ra')
        mixture_table_rows_html = []
        for _, row in top_20_best_pairs.iterrows():
            s1, s2 = row['Solvent 1'], row['Solvent 2']
            pair_data = mixtures_df[(mixtures_df['Solvent 1'] == s1) & (mixtures_df['Solvent 2'] == s2)]
            ra_values_json = json.dumps(list(pair_data.sort_values('Vol Frac 1')['Ra'])); slider_start_index = int(round((row['Vol Frac 1']/0.05)-1))
            mixture_table_rows_html.append(f"""<tr data-ra-values='{ra_values_json}'><td>{s1}</td><td>{s2}</td><td>{row['Vol Frac 1']:.0%} / {1-row['Vol Frac 1']:.0%}</td><td>{row['Ra']:.2f}</td><td class="slider-cell"><input type="range" min="0" max="18" value="{slider_start_index}" class="ra-slider"><span class="slider-display"></span></td></tr>""")
        html_output = f"""<h3 class="table-title">{title}</h3><table class="interactive-table"><thead><tr><th>Solvent 1</th><th>Solvent 2</th><th>Optimized Ratio</th><th>Lowest Ra</th><th>Explore Ratios</th></tr></thead><tbody>{''.join(mixture_table_rows_html)}</tbody></table>"""
        if is_detailed_view: return f"<details><summary>Show optimized systems using all available solvents</summary>{html_output}</details>"
        return html_output

    top_common_mixtures_html = get_optimized_mixtures_html(common_solvents_only_df, "Using Common Solvents")
    top_all_mixtures_html = get_optimized_mixtures_html(solvent_df, "Using All Solvents", is_detailed_view=True)

    target_hsp_json = json.dumps({"d": final_delta_d, "p": final_delta_p, "h": final_delta_h})
    solvents_json = solvent_df[['Solvent', 'δD', 'δP', 'δH']].to_json(orient='records')

    # --- Part 3: HTML Assembly ---
    mol_image = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol_image, size=(300, 300)); buf = BytesIO(); img.save(buf, format="PNG")
    img_base64 = base64.b64encode(buf.getvalue()).decode("utf-8")
    molecule_structure_html = f"""<div class="structure-card"><img src="data:image/png;base64,{img_base64}" alt="Molecule Structure"/></div>"""

    group_contrib_df = pd.DataFrame.from_dict(group_contributions, orient='index').reset_index()
    group_contrib_df.rename(columns={"index":"Functional Group", "delta_d":"ΣD", "delta_p":"ΣP", "delta_h":"ΣH", "count":"Count"}, inplace=True)
    contribution_tables_html = group_contrib_df[['Functional Group', 'Count', 'ΣD', 'ΣP', 'ΣH']].to_html(index=False, classes='small-table')

    calculated_df = pd.DataFrame({"Parameter": ["δD", "δP", "δH"], "Value": [round(v, 2) for v in [final_delta_d, final_delta_p, final_delta_h]]})

    fig = go.Figure()
    df_common = solvent_df[solvent_df['is_common']]; df_other = solvent_df[~solvent_df['is_common']]
    fig.add_trace(go.Scatter3d(x=df_other["δD"], y=df_other["δP"], z=df_other["δH"], mode="markers", marker=dict(size=5, color='gray', opacity=0.5), text=df_other["Solvent"], hoverinfo="text", name="Other Solvents"))
    fig.add_trace(go.Scatter3d(x=df_common["δD"], y=df_common["δP"], z=df_common["δH"], mode="markers", marker=dict(size=6, color='#1f77b4'), text=df_common["Solvent"], hoverinfo="text", name="Common Solvents"))
    fig.add_trace(go.Scatter3d(x=[final_delta_d], y=[final_delta_p], z=[final_delta_h], mode="markers", marker=dict(size=10, color='purple', symbol='diamond'), name="Target Compound", hovertext=[f"Target<br>δD:{final_delta_d:.2f}<br>δP:{final_delta_p:.2f}<br>δH:{final_delta_h:.2f}"], hoverinfo="text"))
    fig.update_layout(title="Hansen Solubility 3D Plot", scene=dict(xaxis_title='δD (Dispersion)', yaxis_title='δP (Polar)', zaxis_title='δH (H-Bonding)'), margin=dict(l=0, r=0, b=0, t=40), legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
    plot_html = pio.to_html(fig, include_plotlyjs=True, full_html=False)

    html_template = f"""
    <html><head><title>Hansen Solubility Analysis</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; margin: 2em; color: #333; }}
        .report-section {{ margin-bottom: 2.5em; }}
        .side-by-side-container {{ display: flex; gap: 30px; align-items: flex-start; flex-wrap: wrap; }}
        .side-by-side-container > div {{ flex: 1; min-width: 400px; }}
        .structure-card {{ text-align: center; max-width: 320px; margin: auto; padding: 1em; border: 1px solid #ccc; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 1em; border: 1px solid #ccc; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; vertical-align: middle; }}
        th {{ background-color: #f2f2f2; }}
        h1, h2 {{ border-bottom: 2px solid #eee; padding-bottom: 5px; color: #444; margin-top: 1.5em; }}
        h3.table-title {{ margin-top: 1em; border-bottom: 1px dashed #ccc; padding-bottom: 5px; }}
        summary {{ cursor: pointer; font-weight: bold; margin-top: 1em; padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd; border-radius: 4px;}}
        details > summary {{ margin-bottom: 1em; }}
        code {{ background-color: #eee; padding: 2px 4px; border-radius: 3px; }}
        .slider-cell {{ display: flex; align-items: center; gap: 10px; min-width:300px; }}
        .ra-slider {{ flex-grow: 1; }}
        .slider-display {{ font-weight: bold; font-family: monospace; min-width: 130px; }}
        .sandbox-controls {{ display:flex; gap:20px; margin-bottom:1em; align-items: center; }}
        #sandbox-result {{ border: 1px solid #007bff; border-radius: 5px; padding: 15px; margin-top: 1em; display: none; }}
    </style>
    </head><body>
    <div class="report-section"><h1>Hansen Solubility Analysis for <code>{smiles}</code></h1></div>
    <div class="report-section"><h2>Target Compound Structure</h2>{molecule_structure_html}</div>
    <div class="report-section"><h2>Final Hansen Parameters</h2><div>{calculated_df.to_html(index=False, classes='small-table')}</div></div>
    <div class="report-section">
        <h2>Single Solvent Matching</h2>
        <div class="side-by-side-container">
            <div><h3 class="table-title">Top 5 Best Common Solvents</h3>{best_matches.to_html(index=False)}<details><summary>Show complete list (all solvents)</summary>{all_matches.to_html(index=False)}</details></div>
            <div><h3 class="table-title">Top 5 Worst Common Solvents</h3>{worst_matches.to_html(index=False)}</div>
        </div>
    </div>
    <div class="report-section">
        <h2>Hansen Solubility 3D Plot</h2>
        <details><summary>Click to show or hide the interactive plot</summary><div>{plot_html}</div></details>
    </div>
    <div class="report-section">
        <h2>Choose two solvents to find the optimal ratio:</h2>
        <div class="sandbox-controls">
            <label>Solvent 1: <select id="solvent1-select"><option value="">-- Select --</option></select></label>
            <label>Solvent 2: <select id="solvent2-select"><option value="">-- Select --</option></select></label>
        </div>
        <div id="sandbox-result"><p><strong>Optimal Mix:</strong> <span id="sandbox-optimal"></span></p><div class="slider-cell"><input type="range" min="0" max="18" value="9" class="ra-slider" id="sandbox-slider"><span class="slider-display" id="sandbox-display"></span></div></div>
    </div>
    <div class="report-section"><h2>Optimized Two-Component Solvent Systems</h2>{top_common_mixtures_html}{top_all_mixtures_html}</div>
    <div class="report-section"><h2>Functional Group Contribution Details</h2><div>{contribution_tables_html}</div></div>
    <script>
    document.addEventListener('DOMContentLoaded', function() {{
        const targetHSP = {target_hsp_json};
        const allSolvents = {solvents_json};
        function updateSliderDisplay(slider) {{
            const row = slider.closest('tr, div#sandbox-result'); if (!row) return;
            const display = row.querySelector('.slider-display');
            const raValues = JSON.parse(row.dataset.raValues);
            const index = parseInt(slider.value, 10);
            const percent1 = (index + 1) * 5;
            const ra = raValues[index];
            display.textContent = `${{percent1}}% / ${{100-percent1}}% | Ra: ${{ra.toFixed(2)}}`;
        }}
        document.querySelectorAll('.ra-slider:not(#sandbox-slider)').forEach(slider => {{
            updateSliderDisplay(slider);
            slider.addEventListener('input', () => updateSliderDisplay(slider));
        }});
        const s1_select = document.getElementById('solvent1-select'), s2_select = document.getElementById('solvent2-select');
        const sandboxResultDiv = document.getElementById('sandbox-result'), sandboxSlider = document.getElementById('sandbox-slider');
        const sandboxOptimalSpan = document.getElementById('sandbox-optimal');
        allSolvents.forEach(s => {{ s1_select.options.add(new Option(s.Solvent, s.Solvent)); s2_select.options.add(new Option(s.Solvent, s.Solvent)); }});
        function runSandboxCalculation() {{
            const s1_name = s1_select.value, s2_name = s2_select.value;
            if (!s1_name || !s2_name || s1_name === s2_name) {{ sandboxResultDiv.style.display = 'none'; return; }}
            const s1 = allSolvents.find(s => s.Solvent === s1_name), s2 = allSolvents.find(s => s.Solvent === s2_name);
            let calculatedRa = [], minRa = Infinity, optimalIndex = -1;
            for (let i = 0; i < 19; i++) {{
                const volFrac1 = (i + 1) * 0.05;
                const mixD = s1['δD'] * volFrac1 + s2['δD'] * (1-volFrac1), mixP = s1['δP'] * volFrac1 + s2['δP'] * (1-volFrac1), mixH = s1['δH'] * volFrac1 + s2['δH'] * (1-volFrac1);
                const ra = Math.sqrt(4 * (targetHSP.d - mixD)**2 + (targetHSP.p - mixP)**2 + (targetHSP.h - mixH)**2);
                calculatedRa.push(ra);
                if (ra < minRa) {{ minRa = ra; optimalIndex = i; }}
            }}
            sandboxResultDiv.dataset.raValues = JSON.stringify(calculatedRa);
            sandboxSlider.value = optimalIndex;
            const optimalPercent = (optimalIndex + 1) * 5;
            sandboxOptimalSpan.textContent = `${{optimalPercent}}% ${{s1_name}} / ${{100-optimalPercent}}% ${{s2_name}} (Lowest Ra: ${{minRa.toFixed(2)}})`;
            updateSliderDisplay(sandboxSlider);
            sandboxResultDiv.style.display = 'block';
        }}
        s1_select.addEventListener('change', runSandboxCalculation);
        s2_select.addEventListener('change', runSandboxCalculation);
        sandboxSlider.addEventListener('input', () => updateSliderDisplay(sandboxSlider));
    }});
    </script>
    </body></html>"""

    print(f"HTML report generated for {smiles}.")
    return html_template

def plot_and_calculate_for_copolymer(smiles_percentage_dictionary):
    """
    Performs a full analysis for a copolymer and generates a highly interactive and
    well-formatted HTML report.
    """
    # This list can be customized as needed.
    common_solvent_list = [
        'acetic acid', 'acetone', 'acetonitrile', 'benzene', '1-butanol', '2-butanol',
        '2-butanone', 't-butyl alcohol', 'carbon tetrachloride', 'chlorobenzene', 'chloroform',
        'cyclohexane', '1,2-dichloroethane', 'diethylene glycol', 'diethyl ether', 'diglyme',
        '1,2-dimethoxyethane', 'dimethylformamide', 'dimethyl sulfoxide', '1,4-dioxane',
        'ethanol', 'ethyl acetate', 'ethylene glycol', 'glycerin', 'heptane',
        'hexamethylphosphoramide', 'hexane', 'methanol', 'methyl t-butyl ether',
        'methylene chloride', 'n-methyl-2-pyrrolidinone', 'nitromethane', 'pentane',
        '1-propanol', '2-propanol', 'pyridine', 'tetrahydrofuran', 'toluene',
        'triethyl amine', 'water', 'o-xylene', 'm-xylene', 'p-xylene'
    ]

    # --- Part 1: Copolymer Calculation ---
    final_delta_d, final_delta_p, final_delta_h, monomer_details = calculate_hansen_for_copolymer(smiles_percentage_dictionary)
    if final_delta_d is None: return
    print(f"Final Copolymer Hansen Parameters (δD, δP, δH): ({final_delta_d:.2f}, {final_delta_p:.2f}, {final_delta_h:.2f})")

    # --- Part 2: Data Preparation for Report ---
    try:
        solvent_df = pd.read_csv('Solvent_values_hansen_solubility.csv')
    except FileNotFoundError:
        print("Error: 'Solvent_values_hansen_solubility.csv' not found. Aborting.")
        return

    solvent_df['is_common'] = solvent_df['Solvent'].str.lower().isin([s.lower() for s in common_solvent_list])

    def hansen_distance(row):
        return np.sqrt(4 * (final_delta_d - row["δD"])**2 + (final_delta_p - row["δP"])**2 + (final_delta_h - row["δH"])**2)
    solvent_df["Ra"] = solvent_df.apply(hansen_distance, axis=1)

    # Filter for single solvent matching tables
    all_matches = solvent_df.sort_values("Ra")[["Solvent", "Ra", "δD", "δP", "δH", "is_common"]]
    common_solvents_only_df = solvent_df[solvent_df['is_common']]
    best_matches = common_solvents_only_df.nsmallest(5, "Ra")[["Solvent", "Ra", "δD", "δP", "δH"]]
    worst_matches = common_solvents_only_df.nlargest(5, "Ra")[["Solvent", "Ra", "δD", "δP", "δH"]]

    # Helper function for two-component optimization
    def get_optimized_mixtures_html(solvents_to_use_df, title, is_detailed_view=False):
        print(f"Performing optimization for: {title}...")
        mixture_results = []
        if len(solvents_to_use_df) < 2:
            return f"<h3>{title}</h3><p>Not enough solvents to create mixtures.</p>"

        for (idx1, idx2) in combinations(solvents_to_use_df.index, 2):
            s1_data, s2_data = solvents_to_use_df.loc[idx1], solvents_to_use_df.loc[idx2]
            for vol_frac1 in np.arange(0.05, 1.0, 0.05):
                mix_dD = (s1_data['δD'] * vol_frac1) + (s2_data['δD'] * (1.0 - vol_frac1))
                mix_dP = (s1_data['δP'] * vol_frac1) + (s2_data['δP'] * (1.0 - vol_frac1))
                mix_dH = (s1_data['δH'] * vol_frac1) + (s2_data['δH'] * (1.0 - vol_frac1))
                ra_mix = np.sqrt(4*(final_delta_d-mix_dD)**2 + (final_delta_p-mix_dP)**2 + (final_delta_h-mix_dH)**2)
                mixture_results.append({"Solvent 1": s1_data['Solvent'], "Solvent 2": s2_data['Solvent'], "Vol Frac 1": vol_frac1, "Ra": ra_mix})

        if not mixture_results:
             return f"<h3>{title}</h3><p>No mixtures could be generated.</p>"

        mixtures_df = pd.DataFrame(mixture_results)
        best_pairs_df = mixtures_df.loc[mixtures_df.groupby(['Solvent 1', 'Solvent 2'])['Ra'].idxmin()]
        top_20_best_pairs = best_pairs_df.nsmallest(20, 'Ra')
        mixture_table_rows_html = []
        for _, row in top_20_best_pairs.iterrows():
            s1, s2 = row['Solvent 1'], row['Solvent 2']
            pair_data = mixtures_df[(mixtures_df['Solvent 1'] == s1) & (mixtures_df['Solvent 2'] == s2)]
            ra_values_json = json.dumps(list(pair_data.sort_values('Vol Frac 1')['Ra']))
            slider_start_index = int(round((row['Vol Frac 1'] / 0.05) - 1))
            mixture_table_rows_html.append(f"""<tr data-ra-values='{ra_values_json}'><td>{s1}</td><td>{s2}</td><td>{row['Vol Frac 1']:.0%} / {1-row['Vol Frac 1']:.0%}</td><td>{row['Ra']:.2f}</td><td class="slider-cell"><input type="range" min="0" max="18" value="{slider_start_index}" class="ra-slider"><span class="slider-display"></span></td></tr>""")

        html_output = f"""<h3 class="table-title">{title}</h3><table class="interactive-table"><thead><tr><th>Solvent 1</th><th>Solvent 2</th><th>Optimized Ratio</th><th>Lowest Ra</th><th>Explore Ratios</th></tr></thead><tbody>{''.join(mixture_table_rows_html)}</tbody></table>"""

        if is_detailed_view:
            return f"<details><summary>Show optimized systems using all available solvents</summary>{html_output}</details>"
        return html_output

    top_common_mixtures_html = get_optimized_mixtures_html(common_solvents_only_df, "Using Common Solvents")
    top_all_mixtures_html = get_optimized_mixtures_html(solvent_df, "Using All Solvents", is_detailed_view=True)

    copolymer_hsp_json = json.dumps({"d": final_delta_d, "p": final_delta_p, "h": final_delta_h})
    solvents_json = solvent_df[['Solvent', 'δD', 'δP', 'δH']].to_json(orient='records')

    # --- Part 3: HTML Assembly ---
    monomer_html_parts, contribution_tables_html = [], ""
    for smiles, details in monomer_details.items():
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol, size=(250, 250)); buf = BytesIO(); img.save(buf, format="PNG")
        img_base64 = base64.b64encode(buf.getvalue()).decode("utf-8")
        monomer_html_parts.append(f"""<div class="monomer-card"><h3>Copolymer section: <code>{smiles}</code></h3><p><b>Mole Fraction: {details['mole_fraction']:.2%}</b></p><img src="data:image/png;base64,{img_base64}" alt="Monomer structure for {smiles}" /></div>""")
        group_contrib_df = pd.DataFrame.from_dict(details['group_contributions'], orient='index').reset_index()
        if not group_contrib_df.empty:
            group_contrib_df.rename(columns={"index": "Functional Group", "delta_d": "ΣD", "delta_p": "ΣP", "delta_h": "ΣH", "count": "Count"}, inplace=True)
            contribution_tables_html += f"<h4>Contributions from <code>{smiles}</code> ({details['mole_fraction']:.2%})</h4>" + group_contrib_df[['Functional Group', 'Count', 'ΣD', 'ΣP', 'ΣH']].to_html(index=False, classes='small-table')
    monomers_section_html = "".join(monomer_html_parts)
    calculated_df = pd.DataFrame({"Parameter": ["Final Copolymer δD", "Final Copolymer δP", "Final Copolymer δH"], "Value": [round(v, 2) for v in [final_delta_d, final_delta_p, final_delta_h]]})

    fig = go.Figure()
    df_common = solvent_df[solvent_df['is_common']]
    df_other = solvent_df[~solvent_df['is_common']]
    fig.add_trace(go.Scatter3d(x=df_other["δD"], y=df_other["δP"], z=df_other["δH"], mode="markers", marker=dict(size=5, color='gray', opacity=0.5), text=df_other["Solvent"], hoverinfo="text", name="Other Solvents"))
    fig.add_trace(go.Scatter3d(x=df_common["δD"], y=df_common["δP"], z=df_common["δH"], mode="markers", marker=dict(size=6, color='#1f77b4'), text=df_common["Solvent"], hoverinfo="text", name="Common Solvents"))
    fig.add_trace(go.Scatter3d(x=[final_delta_d], y=[final_delta_p], z=[final_delta_h], mode="markers", marker=dict(size=10, color='purple', symbol='diamond'), name="Copolymer", hovertext=[f"Copolymer<br>δD:{final_delta_d:.2f}<br>δP:{final_delta_p:.2f}<br>δH:{final_delta_h:.2f}"], hoverinfo="text"))
    fig.update_layout(title="Hansen Solubility 3D Plot", scene=dict(xaxis_title='δD (Dispersion)', yaxis_title='δP (Polar)', zaxis_title='δH (H-Bonding)'), margin=dict(l=0, r=0, b=0, t=40), legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
    plot_html = pio.to_html(fig, include_plotlyjs=True, full_html=False)

    html_template = f"""
    <html><head><title>Hansen Solubility Analysis</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; margin: 2em; color: #333; }}
        .report-section {{ margin-bottom: 2.5em; }}
        .side-by-side-container {{ display: flex; gap: 30px; align-items: flex-start; flex-wrap: wrap; }}
        .side-by-side-container > div {{ flex: 1; min-width: 400px; }}
        .monomer-section {{ display: flex; flex-wrap: wrap; gap: 20px; width: 100%; }}
        .monomer-card {{ flex: 1; min-width: 280px; border: 1px solid #ccc; border-radius: 8px; padding: 1em; box-shadow: 0 2px 5px rgba(0,0,0,0.1); text-align: center; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 1em; border: 1px solid #ccc; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; vertical-align: middle; }}
        th {{ background-color: #f2f2f2; }}
        h1, h2 {{ border-bottom: 2px solid #eee; padding-bottom: 5px; color: #444; margin-top: 1.5em; }}
        h3.table-title {{ margin-top: 1em; border-bottom: 1px dashed #ccc; padding-bottom: 5px; }}
        summary {{ cursor: pointer; font-weight: bold; margin-top: 1em; padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd; border-radius: 4px;}}
        details > summary {{ margin-bottom: 1em; }}
        code {{ background-color: #eee; padding: 2px 4px; border-radius: 3px; }}
        .slider-cell {{ display: flex; align-items: center; gap: 10px; min-width:300px; }}
        .ra-slider {{ flex-grow: 1; }}
        .slider-display {{ font-weight: bold; font-family: monospace; min-width: 130px; }}
        .sandbox-controls {{ display:flex; gap:20px; margin-bottom:1em; align-items: center; }}
        #sandbox-result {{ border: 1px solid #007bff; border-radius: 5px; padding: 15px; margin-top: 1em; display: none; }}
    </style>
    </head><body>
    <div class="report-section"><h1>Hansen Solubility Analysis for Copolymer</h1></div>

    <div class="report-section">
        <h2>Copolymer Composition</h2>
        <div class="monomer-section">{monomers_section_html}</div>
    </div>

    <div class="report-section">
        <h2>Final Copolymer Results</h2>
        <div>{calculated_df.to_html(index=False, classes='small-table')}</div>
    </div>

    <div class="report-section">
        <h2>Single Solvent Matching</h2>
        <div class="side-by-side-container">
            <div><h3 class="table-title">Top 5 Best Common Solvents</h3>{best_matches.to_html(index=False)}<details><summary>Show complete list (all solvents)</summary>{all_matches.to_html(index=False)}</details></div>
            <div><h3 class="table-title">Top 5 Worst Common Solvents</h3>{worst_matches.to_html(index=False)}</div>
        </div>
    </div>

    <div class="report-section">
        <h2>Hansen Solubility 3D Plot</h2>
        <details>
            <summary>Click to show or hide the interactive plot</summary>
            <div>{plot_html}</div>
        </details>
    </div>

    <div class="report-section">
        <h2>Choose two solvents to optimise their ratio for maximum solubility:</h2>
        <div class="sandbox-controls">
            <label>Solvent 1: <select id="solvent1-select"><option value="">-- Select --</option></select></label>
            <label>Solvent 2: <select id="solvent2-select"><option value="">-- Select --</option></select></label>
        </div>
        <div id="sandbox-result">
            <p><strong>Optimal Mix:</strong> <span id="sandbox-optimal"></span></p>
            <div class="slider-cell"><input type="range" min="0" max="18" value="9" class="ra-slider" id="sandbox-slider"><span class="slider-display" id="sandbox-display"></span></div>
        </div>
    </div>

    <div class="report-section">
        <h2>Optimized Two-Component Solvent Systems</h2>
        {top_common_mixtures_html}
        {top_all_mixtures_html}
    </div>

    <div class="report-section">
        <h2>Functional Group Contribution Details</h2>
        <div>{contribution_tables_html}</div>
    </div>

    <script>
    document.addEventListener('DOMContentLoaded', function() {{
        const copolymerHSP = {copolymer_hsp_json};
        const allSolvents = {solvents_json};
        function updateSliderDisplay(slider) {{
            const row = slider.closest('tr, div#sandbox-result'); if (!row) return;
            const display = row.querySelector('.slider-display');
            const raValues = JSON.parse(row.dataset.raValues);
            const index = parseInt(slider.value, 10);
            const percent1 = (index + 1) * 5;
            const ra = raValues[index];
            display.textContent = `${{percent1}}% / ${{100-percent1}}% | Ra: ${{ra.toFixed(2)}}`;
        }}
        document.querySelectorAll('.ra-slider:not(#sandbox-slider)').forEach(slider => {{
            updateSliderDisplay(slider);
            slider.addEventListener('input', () => updateSliderDisplay(slider));
        }});
        const s1_select = document.getElementById('solvent1-select'), s2_select = document.getElementById('solvent2-select');
        const sandboxResultDiv = document.getElementById('sandbox-result'), sandboxSlider = document.getElementById('sandbox-slider');
        const sandboxOptimalSpan = document.getElementById('sandbox-optimal');
        allSolvents.forEach(s => {{ s1_select.options.add(new Option(s.Solvent, s.Solvent)); s2_select.options.add(new Option(s.Solvent, s.Solvent)); }});
        function runSandboxCalculation() {{
            const s1_name = s1_select.value, s2_name = s2_select.value;
            if (!s1_name || !s2_name || s1_name === s2_name) {{ sandboxResultDiv.style.display = 'none'; return; }}
            const s1 = allSolvents.find(s => s.Solvent === s1_name), s2 = allSolvents.find(s => s.Solvent === s2_name);
            let calculatedRa = [], minRa = Infinity, optimalIndex = -1;
            for (let i = 0; i < 19; i++) {{
                const volFrac1 = (i + 1) * 0.05;
                const mixD = s1['δD'] * volFrac1 + s2['δD'] * (1-volFrac1), mixP = s1['δP'] * volFrac1 + s2['δP'] * (1-volFrac1), mixH = s1['δH'] * volFrac1 + s2['δH'] * (1-volFrac1);
                const ra = Math.sqrt(4 * (copolymerHSP.d - mixD)**2 + (copolymerHSP.p - mixP)**2 + (copolymerHSP.h - mixH)**2);
                calculatedRa.push(ra);
                if (ra < minRa) {{ minRa = ra; optimalIndex = i; }}
            }}
            sandboxResultDiv.dataset.raValues = JSON.stringify(calculatedRa);
            sandboxSlider.value = optimalIndex;
            const optimalPercent = (optimalIndex + 1) * 5;
            sandboxOptimalSpan.textContent = `${{optimalPercent}}% ${{s1_name}} / ${{100-optimalPercent}}% ${{s2_name}} (Lowest Ra: ${{minRa.toFixed(2)}})`;
            updateSliderDisplay(sandboxSlider);
            sandboxResultDiv.style.display = 'block';
        }}
        s1_select.addEventListener('change', runSandboxCalculation);
        s2_select.addEventListener('change', runSandboxCalculation);
        sandboxSlider.addEventListener('input', () => updateSliderDisplay(sandboxSlider));
    }});
    </script>
    </body></html>"""

    print(f"HTML report generated for copolymer.")
    return html_template