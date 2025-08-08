from flask import Flask, request, render_template_string, jsonify
import traceback
import json

# The import statements now include the new show_structure function
from solubility_prediction import plot_and_calculate_values, show_structure, plot_and_calculate_for_copolymer

# Initialize the Flask app
app = Flask(__name__)


# --- HTML Template for the Input Form ---
# This is the homepage the user will see.
HOME_PAGE_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Hansen Solubility Calculator for Polymers</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;700&display=swap" rel="stylesheet">

    <style>
        :root {
            --bg-color: #f8f9fa;
            --card-bg: #ffffff;
            --text-color: #343a40;
            --heading-color: #212529;
            --primary-color: #007bff;
            --primary-hover: #0056b3;
            --secondary-color: #6c757d;
            --secondary-hover: #5a6268;
            --border-color: #dee2e6;
            --shadow: 0 4px 6px rgba(0, 0, 0, 0.05);
        }
        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            margin: 0;
            padding: 2em;
            background-color: var(--bg-color);
            color: var(--text-color);
            line-height: 1.6;
        }
        .container {
            max-width: 900px;
            margin: 0 auto;
        }
        header {
            text-align: center;
            margin-bottom: 3em;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 1.5em;
        }
        header h1 {
            font-size: 2.5em;
            color: var(--heading-color);
            margin: 0;
        }
        .card {
            background-color: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 2em;
            margin-bottom: 2.5em;
            box-shadow: var(--shadow);
        }
        h1, h2, h3 {
            color: var(--heading-color);
            font-weight: 700;
        }
        h2 {
            margin-top: 0;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 0.5em;
            margin-bottom: 1em;
        }
        input[type="text"], input[type="number"] {
            font-size: 1.1em;
            padding: 12px 15px;
            width: 100%;
            max-width: 600px;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            box-sizing: border-box;
            transition: border-color 0.2s, box-shadow 0.2s;
        }
        input[type="text"]:focus, input[type="number"]:focus {
            outline: none;
            border-color: var(--primary-color);
            box-shadow: 0 0 0 3px rgba(0, 123, 255, 0.25);
        }
        button {
            font-size: 1.1em;
            font-weight: 500;
            padding: 12px 20px;
            border: none;
            border-radius: 6px;
            color: white;
            cursor: pointer;
            transition: background-color 0.2s;
        }
        button.primary {
            background-color: var(--primary-color);
        }
        button.primary:hover {
            background-color: var(--primary-hover);
        }
        button.secondary {
            background-color: var(--secondary-color);
        }
        button.secondary:hover {
            background-color: var(--secondary-hover);
        }
        button.danger {
            background-color: #dc3545;
        }
        button.danger:hover {
            background-color: #c82333;
        }
        button[type="submit"] {
            margin-left: 10px;
        }
        .form-group {
            display: flex;
            align-items: center;
        }
        .action-buttons {
            display: flex;
            justify-content: center;
            gap: 1em;
            margin-top: 1.5em;
        }
        .structure-result {
            margin-top: 2em;
            text-align: center;
            background-color: #fdfdff;
            padding: 1em;
            border-radius: 8px;
            border: 1px dashed var(--border-color);
            min-height: 50px; /* Ensures space is reserved */
        }
        .structure-result code {
            background-color: #e9ecef;
            padding: 2px 5px;
            border-radius: 4px;
        }
        ol, ul {
            padding-left: 20px;
        }
        li {
            margin-bottom: 1em;
        }
        footer {
            text-align: center;
            margin-top: 4em;
            padding-top: 2em;
            border-top: 1px solid var(--border-color);
            color: #6c757d;
        }
        footer ul {
            list-style: none;
            padding: 0;
        }
        footer li {
            margin-bottom: 0.5em;
        }
        details > summary {
            cursor: pointer;
            font-weight: 700;
            padding: 0.5em;
            border-radius: 6px;
            transition: background-color 0.2s;
        }
        details > summary:hover {
            background-color: #e9ecef;
        }
        #copolymer-table {
            width: 100%;
            margin-top: 1em;
            border-collapse: collapse;
        }
        #copolymer-table th, #copolymer-table td {
            border: 1px solid var(--border-color);
            padding: 12px;
            text-align: left;
            vertical-align: middle;
        }
        #copolymer-table th {
            background-color: #f8f9fa;
        }
        #copolymer-table img {
            max-width: 150px;
            display: block;
            margin: auto;
        }
        #copolymer-table .percent-input {
            width: 80px;
        }
        #copolymer-total-row {
            font-weight: bold;
        }
        .tab-buttons {
            display: flex;
            border-bottom: 1px solid var(--border-color);
            margin-bottom: 1.5em;
        }
        .tab-button {
            padding: 10px 20px;
            cursor: pointer;
            border: none;
            background-color: transparent;
            font-size: 1.1em;
            font-weight: 500;
            color: var(--secondary-color);
            border-bottom: 3px solid transparent;
            margin-bottom: -1px;
        }
        .tab-button.active {
            color: var(--primary-color);
            border-bottom: 3px solid var(--primary-color);
        }
        .tab-content {
            display: none;
        }
        .tab-content.active {
            display: block;
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Hansen Solubility Calculator for Polymers</h1>
        </header>

        <main>
            <div class="card">
                <details>
                    <summary><h2 style="display: inline-block; margin-bottom: 0; border-bottom: none; padding-bottom: 0;">How to use this website</h2></summary>
                    <div style="padding-top: 1em;">
                        <ol>
                            <li>Select either <strong>Homopolymer</strong> or <strong>Copolymer</strong> analysis.</li>
                            <li><strong>For Homopolymers:</strong> Draw your structure or paste a SMILES string and click "Run Analysis".</li>
                            <li><strong>For Copolymers:</strong> Draw each monomer, click "Add to Copolymer", set the mole percentages in the builder table, and then click "Run Copolymer Analysis".</li>
                            <li>You can check any structure (drawn or pasted) at any time using the "Check Structure" buttons.</li>
                            <li><strong>View the Report:</strong> A new page will load containing the complete, interactive solubility report for your molecule.</li>
                        <li><strong>Thoroughly check the results:</strong>
                            <p>This is a crudely made program, where I attempt to split compounds into their functional groups based on the SMILES you provide, and therefore perform Hansen Solubility Analysis based on these functional groups. The more complicated the structure, the more likely the program will make errors in splitting the structure into it's functional groups. So you can check this, I have included this in a table at the end of the results page. If this is incorrect, missing parts, or has incorrectly assigned functional groups, the calculations will be wrong. These calculations are also not suitable for small molecules, as the calculations for these are significantly more complicated. Please only use this tool for polymers, and take results with a grain of salt. For full information on how these calculations have been performed, please see the explainer at the bottom of this page.</p>
                        </li>
                        </ol>
                    </div>
                </details>
            </div>

            <!-- Main Analysis Section with Tabs -->
            <div class="card">
                <div class="tab-buttons">
                    <button class="tab-button active" onclick="openTab(event, 'homopolymer-tab')">Homopolymer Analysis</button>
                    <button class="tab-button" onclick="openTab(event, 'copolymer-tab')">Copolymer Analysis</button>
                </div>

                <!-- Homopolymer Tab -->
                <div id="homopolymer-tab" class="tab-content active">
                    <h3>Draw Polymer's Repeating Unit</h3>
                    <p>Use the editor to draw the core structure of your polymer's repeating unit. To define how it connects into a chain, you must mark the two connection points using 'R' atoms.</p>
<ol style="text-align: left; max-width: 600px; margin: 1em auto;">
    <li>Draw the single repeating unit of your polymer.</li>
    <li>At the two points where the unit would bond to its neighbors, add an 'R' atom by hovering over the atom and pressing R on your keyboard.</li>
    <li>These 'R' atoms will be automatically converted to the <code>[*]</code> character required for a valid polymer SMILES string.</li>
</ol>
                    <div id="ocl-editor-container-homo" style="width: 100%; height: 450px; border: 1px solid var(--border-color); border-radius: 6px; position: relative; overflow: hidden;"></div>
                    <div class="action-buttons">
                        <button type="button" id="check_drawn_homo_btn" class="secondary">Check Drawn Structure</button>
                        <form id="analyze_drawn_homo_form" action="/analyze" method="post" style="display: inline;">
                            <input type="hidden" name="smiles_analyze">
                            <button type="submit" class="primary">Run Analysis on Drawn Structure</button>
                        </form>
                    </div>
                    <div class="structure-result" id="drawn_homo_result_container" style="display: none;"></div>
                    <hr style="margin: 2em 0;">
                    <h3>Or, Input SMILES String</h3>
                    <form action="/analyze" method="post" class="form-group">
                        <input type="text" name="smiles_analyze" placeholder="e.g., [*]CC([*])C1=CC=CC=C1 for poly(styrene)" required>
                        <button type="submit" class="primary">Run Analysis</button>
                    </form>
                </div>

                <!-- Copolymer Tab -->
                <div id="copolymer-tab" class="tab-content">
                    <h3>Draw and Add Copolymer Units</h3>
                    <p>Use the editor to draw the core structure of your polymer's repeating unit. To define how it connects into a chain, you must mark the two connection points using 'R' atoms.</p>
<ol style="text-align: left; max-width: 600px; margin: 1em auto;">
    <li>Draw the single repeating unit of your polymer.</li>
    <li>At the two points where the unit would bond to its neighbors, add an 'R' atom by hovering over the atom and pressing R on your keyboard.</li>
    <li>These 'R' atoms will be automatically converted to the <code>[*]</code> character required for a valid polymer SMILES string.</li>
</ol>
                    <div id="ocl-editor-container-copo" style="width: 100%; height: 450px; border: 1px solid var(--border-color); border-radius: 6px; position: relative; overflow: hidden;"></div>
                    <div class="action-buttons">
                         <button type="button" id="check_drawn_copo_btn" class="secondary">Check Drawn Structure</button>
                         <button type="button" id="add_to_copolymer_btn" class="primary">Add to Copolymer Builder</button>
                    </div>
                     <div class="structure-result" id="drawn_copo_result_container" style="display: none;"></div>
                    <hr style="margin: 2em 0;">
                    <h3>2. Build Copolymer</h3>
                    <p>Define the mole percentage for each section, ensuring they sum to 100%.</p>
                    <table id="copolymer-table">
                        <thead>
                            <tr>
                                <th>Structure</th>
                                <th>SMILES</th>
                                <th>Mole %</th>
                                <th>Action</th>
                            </tr>
                        </thead>
                        <tbody id="copolymer-table-body"></tbody>
                        <tfoot>
                            <tr id="copolymer-total-row">
                                <td colspan="2" style="text-align: right;"><strong>Total:</strong></td>
                                <td id="copolymer-total-percent">0%</td>
                                <td></td>
                            </tr>
                        </tfoot>
                    </table>
                    <div class="action-buttons">
                         <form id="analyze_copolymer_form" action="/analyze_copolymer" method="post">
                            <input type="hidden" name="copolymer_data">
                            <button type="submit" class="primary">Run Copolymer Analysis</button>
                        </form>
                    </div>
                </div>
            </div>

            <div class="card">
                <details>
                    <summary><strong>Click to see how the calculations are performed</strong></summary>
                    <div style="line-height: 1.6; padding-top: 1em;">
                        <p>The Hansen Solubility Parameters (HSP) for the input polymer are estimated using a group-contribution method. This approach breaks the polymer's repeating unit down into its fundamental functional groups and sums their individual contributions to determine the final properties.</p>
                        <h4>Step 1: Functional Group Identification</h4>
                        <p>The provided polymer SMILES string is parsed using the RDKit library to identify and count the first-order functional groups present in the repeating unit. The accuracy of this step is critical, and for highly complex or unusual structures, it's important to check the "Functional Group Contribution Details" table in the final report to ensure the groups have been identified correctly.</p>
                        <h4>Step 2: Summation of Group Contributions</h4>
                        <p>Each functional group (i) has known contribution values for the Dispersion (D<sub>i</sub>), Polar (P<sub>i</sub>), and Hydrogen Bonding (H<sub>i</sub>) components. These values are sourced from the work of Stefanis and Panayiotou (2008). The total contribution for each parameter is calculated by summing the contributions of all identified groups, multiplied by their count (N<sub>i</sub>):</p>
                        <ul>
                            <li>&delta;<sub>D</sub> = &Sigma; N<sub>i</sub>D<sub>i</sub></li>
                            <li>&delta;<sub>P</sub> = &Sigma; N<sub>i</sub>P<sub>i</sub></li>
                            <li>&delta;<sub>H</sub> = &Sigma; N<sub>i</sub>H<sub>i</sub></li>
                        </ul>
                        <h4>Step 3: Final Calibration</h4>
                        <p>After the initial summation, a set of specific calibration constants are added to each parameter to yield the final estimated HSP values for the polymer. These values are also sourced from the first order contribution calculations from Stefanis and Panayiotou (2008).</p>
                        <h4>Step 4: Solvent Matching (Hansen Distance)</h4>
                        <p>The final predicted HSP of the polymer (&delta;<sub>D,p</sub>, &delta;<sub>P,p</sub>, &delta;<sub>H,p</sub>) is compared against a database of known solvents (&delta;<sub>D,s</sub>, &delta;<sub>P,s</sub>, &delta;<sub>H,s</sub>). The compatibility is determined by calculating the Hansen Distance (R<sub>a</sub>), where a smaller distance signifies a better match. The formula used is:</p>
                        <p><code>Ra = sqrt(4*(&delta;D<sub>polymer</sub> - &delta;D<sub>solvent</sub>)<sup>2</sup> + (&delta;P<sub>polymer</sub> - &delta;P<sub>solvent</sub>)<sup>2</sup> + (&delta;H<sub>polymer</sub> - &delta;H<sub>solvent</sub>)<sup>2</sup>)</code></p>
                        <p>The factor of 4 applied to the dispersion term is a standard convention to better represent the spherical nature of the solubility volume. Solvent data has been gathered from Hansen's website, as well as the Stenutz website cited below.</p>
                        <h4>Step 5: Two-Component Solvent Blends</h4>
                        <p>For solvent mixtures, the blended HSP is calculated as a volume-weighted average of the two individual solvents before the Hansen Distance (R<sub>a</sub>) is calculated against the polymer.</p>
                        <hr>
                        <p><strong>Disclaimer:</strong> This group-contribution method is designed for estimating the properties of polymers, and can be incorrect. As noted by Hansen and others, it is not suitable for small molecules where molar volume plays a more significant role in the calculations. Please use this tool accordingly.</p>
                    </div>
                </details>
            </div>
        </main>

        <footer>
            <h3>Sources</h3>
            <ul>
                <li>Hansen solubility parameters | hansen solubility parameters. <a href="https://www.hansen-solubility.com/">https://www.hansen-solubility.com/</a> (accessed 2025-08-06).</li>
                <li>Hansen solubility parameters. <a href="https://www.stenutz.eu/chem/hansen.php/">https://www.stenutz.eu/chem/hansen.php/</a> (accessed 2025-08-06).</li>
                <li>Stefanis, E.; Panayiotou, C. Prediction of Hansen Solubility Parameters with a New Group-Contribution Method. Int J Thermophys 2008, 29 (2), 568–585. <a href="https://doi.org/10.1007/s10765-008-0415-z">https://doi.org/10.1007/s10765-008-0415-z/</a>.</li>
            </ul>
        </footer>
        <footer>
            <h3>Feedback</h3>
            <ul>
                <li>Any feedback is welcome. If you have had any issues with calculations, have suggested improvements, or have a compound that cannot be split into the correct functional groups, please contact me by email.</li>
                <li><a href="mailto:theo.ellingsen@utas.edu.au">theo.ellingsen@utas.edu.au</a></li>
                <li></li>
            </ul>
        </footer>
    </div>

    <script src="https://www.lactame.com/lib/openchemlib/HEAD/openchemlib-full.pretty.js"></script>

    <script>
        let editorHomo, editorCopo;

        function openTab(evt, tabName) {
            var i, tabcontent, tabbuttons;
            tabcontent = document.getElementsByClassName("tab-content");
            for (i = 0; i < tabcontent.length; i++) {
                tabcontent[i].style.display = "none";
            }
            tabbuttons = document.getElementsByClassName("tab-button");
            for (i = 0; i < tabbuttons.length; i++) {
                tabbuttons[i].className = tabbuttons[i].className.replace(" active", "");
            }
            document.getElementById(tabName).style.display = "block";
            evt.currentTarget.className += " active";

            // Initialize editors if they don't exist yet
            if (tabName === 'homopolymer-tab' && !editorHomo && typeof OCL !== 'undefined') {
                editorHomo = OCL.StructureEditor.createEditor('ocl-editor-container-homo');
                editorHomo.setSmiles('[*]CC([*])C1=CC=CC=C1');
            }
            if (tabName === 'copolymer-tab' && !editorCopo && typeof OCL !== 'undefined') {
                editorCopo = OCL.StructureEditor.createEditor('ocl-editor-container-copo');
                editorCopo.setSmiles('[*]CC([*])C(O)=O');
            }
        }

        document.addEventListener('DOMContentLoaded', function () {
            // Open the default tab and initialize its editor
            document.querySelector('.tab-button.active').click();

            function checkStructure(smilesString, targetContainerId) {
                const targetContainer = document.getElementById(targetContainerId);
                if (!smilesString) {
                    targetContainer.innerHTML = "<p style='color: orange;'>Please provide a SMILES string.</p>";
                    targetContainer.style.display = 'block';
                    return;
                }
                targetContainer.innerHTML = "<p>Loading...</p>";
                targetContainer.style.display = 'block';

                fetch('/check_structure_ajax', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles_check: smilesString }),
                })
                .then(response => response.json())
                .then(data => {
                    if (data.success) {
                        targetContainer.innerHTML = `
                            <h2>Structure Preview</h2>
                            <img src="data:image/png;base64,${data.image_data}" alt="Molecule structure" />
                            <p><code>${data.smiles}</code></p>`;
                    } else {
                        targetContainer.innerHTML = `<p style='color: red;'>${data.error}</p>`;
                    }
                });
            }

            // --- Homopolymer Logic ---
            document.getElementById('check_drawn_homo_btn').addEventListener('click', () => {
                if (editorHomo) {
                    let smiles = editorHomo.getSmiles().replace(/\[R\d*\]/g, '[*]');
                    checkStructure(smiles, 'drawn_homo_result_container');
                }
            });

            document.getElementById('analyze_drawn_homo_form').addEventListener('submit', (event) => {
                if (editorHomo) {
                    let smiles = editorHomo.getSmiles().replace(/\[R\d*\]/g, '[*]');
                    event.target.querySelector('input[name="smiles_analyze"]').value = smiles;
                } else {
                    event.preventDefault();
                    alert("Editor not loaded.");
                }
            });

            // --- Copolymer Logic ---
            const addCopolymerBtn = document.getElementById('add_to_copolymer_btn');
            const copolymerTableBody = document.getElementById('copolymer-table-body');
            const totalPercentCell = document.getElementById('copolymer-total-percent');
            const analyzeCopolymerForm = document.getElementById('analyze_copolymer_form');

            document.getElementById('check_drawn_copo_btn').addEventListener('click', () => {
                 if (editorCopo) {
                    let smiles = editorCopo.getSmiles().replace(/\[R\d*\]/g, '[*]');
                    checkStructure(smiles, 'drawn_copo_result_container');
                }
            });

            function updateTotalPercentage() {
                let total = 0;
                copolymerTableBody.querySelectorAll('.percent-input').forEach(input => {
                    total += parseFloat(input.value) || 0;
                });
                totalPercentCell.textContent = `${total.toFixed(0)}%`;
                totalPercentCell.style.color = Math.abs(total - 100) < 0.1 ? 'green' : 'red';
            }

            addCopolymerBtn.addEventListener('click', function() {
                if (!editorCopo) return;
                let smiles = editorCopo.getSmiles().replace(/\[R\d*\]/g, '[*]');
                if (smiles) {
                    if (copolymerTableBody.querySelector(`tr[data-smiles="${smiles}"]`)) {
                        alert("This monomer has already been added.");
                        return;
                    }
                    fetch('/check_structure_ajax', {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ smiles_check: smiles }),
                    })
                    .then(response => response.json())
                    .then(data => {
                        if (data.success) {
                            const newRow = document.createElement('tr');
                            newRow.dataset.smiles = smiles;
                            newRow.innerHTML = `
                                <td><img src="data:image/png;base64,${data.image_data}" alt="Monomer structure"></td>
                                <td><code>${smiles}</code></td>
                                <td><input type="number" class="percent-input" min="0" max="100" step="1" value="0"></td>
                                <td><button type="button" class="danger remove-btn">Remove</button></td>
                            `;
                            copolymerTableBody.appendChild(newRow);
                            newRow.querySelector('.percent-input').addEventListener('input', updateTotalPercentage);
                            newRow.querySelector('.remove-btn').addEventListener('click', () => {
                                newRow.remove();
                                updateTotalPercentage();
                            });
                            updateTotalPercentage();
                        } else {
                            alert(`Error adding monomer: ${data.error}`);
                        }
                    });
                }
            });

            analyzeCopolymerForm.addEventListener('submit', function(event) {
                event.preventDefault();
                let total = 0;
                const copolymerData = {};
                const rows = copolymerTableBody.querySelectorAll('tr');
                if (rows.length === 0) {
                    alert("Please add at least one monomer.");
                    return;
                }
                rows.forEach(row => {
                    const smiles = row.dataset.smiles;
                    const percentage = parseFloat(row.querySelector('.percent-input').value) || 0;
                    if (percentage > 0) {
                        copolymerData[smiles] = percentage;
                        total += percentage;
                    }
                });
                if (Math.abs(total - 100) > 0.1) {
                    alert(`Total mole percentage must be 100%, but it is ${total.toFixed(0)}%.`);
                    return;
                }
                analyzeCopolymerForm.querySelector('input[name="copolymer_data"]').value = JSON.stringify(copolymerData);
                analyzeCopolymerForm.submit();
            });
        });
    </script>
    </body>
</html>
"""


# --- Define the Routes (Webpages) ---
@app.route('/', methods=['GET', 'POST'])
def home():
    return render_template_string(HOME_PAGE_HTML)

@app.route('/check_structure_ajax', methods=['POST'])
def check_structure_ajax():
    data = request.get_json()
    smiles_string = data.get('smiles_check', '')
    if not smiles_string:
        return jsonify({'success': False, 'error': 'No SMILES string provided.'})
    try:
        img_base64 = show_structure(smiles_string, False)
        if img_base64:
            return jsonify({'success': True, 'image_data': img_base64, 'smiles': smiles_string})
        else:
            return jsonify({'success': False, 'error': 'Invalid SMILES string provided. Please try again.'})
    except Exception as e:
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'An server error occurred while generating the structure.'})

@app.route('/analyze', methods=['POST'])
def analyze():
    smiles_string = request.form.get('smiles_analyze', '')
    if not smiles_string:
        return "Error: No SMILES string provided for analysis."
    try:
        report_html = plot_and_calculate_values(smiles_string)
        return report_html
    except Exception:
        error_details = traceback.format_exc()
        return f"<h1>An Error Occurred</h1><p>The server encountered an error while processing your request. Here are the details:</p><pre>{error_details}</pre>"

@app.route('/analyze_copolymer', methods=['POST'])
def analyze_copolymer():
    copolymer_data_str = request.form.get('copolymer_data', '{}')
    try:
        # Convert percentages from strings to floats
        smiles_percentage_dictionary = {k: float(v) for k, v in json.loads(copolymer_data_str).items()}
        if not smiles_percentage_dictionary:
             return "Error: No copolymer data provided for analysis."
        report_html = plot_and_calculate_for_copolymer(smiles_percentage_dictionary)
        return report_html
    except Exception:
        error_details = traceback.format_exc()
        return f"<h1>An Error Occurred During Copolymer Analysis</h1><p>The server encountered an error. Here are the details:</p><pre>{error_details}</pre>"


# ADDED: Block to run the app locally for testing
#if __name__ == '__main__':
#    app.run(debug=True, port=50