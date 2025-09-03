"""
Genomics Surveillance Web App (Flask + Biopython)
- Upload or paste FASTA sequences
- Detect substitutions, indels
- Compute a simple heuristic "risk score"

Run:
    python app.py
Then open http://127.0.0.1:5000 in Chrome
"""

from flask import Flask, request, render_template_string, jsonify
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from io import StringIO

app = Flask(__name__)

# ---------------- HTML Template ---------------- #
INDEX_HTML = """
<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Genomics Surveillance — Simple Demo</title>
  <style>
    body{font-family: Arial, sans-serif; max-width:900px;margin:auto;padding:20px}
    textarea{width:100%;height:120px}
    table{width:100%;border-collapse:collapse;margin-top:12px}
    th,td{border:1px solid #ddd;padding:6px;text-align:left}
    th{background:#f4f4f4}
    .mut{color:#b22222}
  </style>
</head>
<body>
  <h1>Genomics Surveillance — Demo</h1>
  <p>Upload a <strong>reference</strong> FASTA and one <strong>sample</strong> FASTA (or paste sequences).</p>

  <form id="seqForm" method="post" action="/analyze" enctype="multipart/form-data">
    <label>Reference (FASTA file or paste):</label><br>
    <input type="file" name="ref_file"><br><small>or paste below:</small>
    <textarea name="ref_text" placeholder=">ref\nATGC..."></textarea>

    <label>Sample (FASTA file or paste):</label><br>
    <input type="file" name="sample_file"><br><small>or paste below:</small>
    <textarea name="sample_text" placeholder=">sample\nATGC..."></textarea>

    <p><button type="submit">Analyze</button></p>
  </form>

  <div id="result"></div>

<script>
const form = document.getElementById('seqForm');
form.onsubmit = async (e)=>{
  e.preventDefault();
  const fd = new FormData(form);
  const res = await fetch('/analyze',{method:'POST',body:fd});
  const data = await res.json();
  const out = document.getElementById('result');
  if(data.error){ out.innerHTML = '<pre style="color:red">'+data.error+'</pre>'; return }
  let html = `<h2>Summary</h2>
  <p>Total variants: <strong>${data.summary.total_variants}</strong></p>
  <p>Substitutions: ${data.summary.substitutions} | Indels: ${data.summary.indels}</p>
  <p>Risk score (0-100): <strong>${data.summary.risk_score}</strong></p>`;

  html += '<h3>Variants</h3>';
  html += '<table><thead><tr><th>Pos</th><th>Ref</th><th>Alt</th><th>Type</th></tr></thead><tbody>';
  data.variants.forEach(v=>{
    html += `<tr><td>${v.pos}</td><td>${v.ref}</td><td class="mut">${v.alt}</td><td>${v.type}</td></tr>`;
  });
  html += '</tbody></table>';

  out.innerHTML = html;
}
</script>
</body>
</html>
"""

# ---------------- Helper Functions ---------------- #

def parse_fasta_from_upload(file_storage, text_field):
    """
    Reads a FASTA sequence either from an uploaded file or from pasted text.
    Handles encoding issues gracefully.
    """
    content = ""
    if file_storage and file_storage.filename:
        try:
            content = file_storage.read().decode("utf-8")
        except UnicodeDecodeError:
            # fallback for Windows-1252 / Latin-1 encodings
            content = file_storage.read().decode("latin-1", errors="ignore")
    else:
        content = text_field or ""

    content = content.strip()
    if not content:
        return None

    handle = StringIO(content)
    try:
        recs = list(SeqIO.parse(handle, "fasta"))
        if len(recs) == 0:
            seq = "".join(content.split())
            return Seq(seq.upper())
        return recs[0].seq
    except Exception:
        seq = "".join(content.split())
        return Seq(seq.upper())


def align_and_compare(ref_seq, sample_seq):
    ref = str(ref_seq)
    sample = str(sample_seq)
    if len(ref) == len(sample):
        return compare_same_length(ref, sample)
    aln = pairwise2.align.globalxx(ref, sample, one_alignment_only=True)[0]
    aligned_ref, aligned_sample = aln.seqA, aln.seqB
    return compare_same_length(aligned_ref, aligned_sample, aligned=True)


def compare_same_length(ref, sample, aligned=False):
    variants = []
    substitutions = 0
    indels = 0
    pos = 0
    for r, s in zip(ref, sample):
        if r == s:
            if r != "-":
                pos += 1
            continue
        if r == "-" or s == "-":
            indels += 1
            var_type = "indel"
        else:
            substitutions += 1
            var_type = "substitution"
        display_pos = pos + 1 if r != "-" else pos + 1
        variants.append({"pos": display_pos, "ref": r, "alt": s, "type": var_type})
        if r != "-":
            pos += 1
    return variants, substitutions, indels


def compute_risk_score(variants, substitutions, indels, ref_len):
    if ref_len == 0:
        return 0
    base = (substitutions + 2 * indels) / max(1, ref_len)
    positions = [v["pos"] for v in variants]
    cluster_factor = 1.0
    if len(positions) >= 2:
        positions.sort()
        max_within = 0
        for i in range(len(positions)):
            j = i
            while j < len(positions) and positions[j] - positions[i] <= 50:
                j += 1
            max_within = max(max_within, j - i)
        cluster_factor = 1.0 + (max_within / max(1, len(positions)))
    raw = base * cluster_factor * 200
    score = max(0, min(100, int(raw)))
    return score

# ---------------- Flask Routes ---------------- #

@app.route("/")
def index():
    return render_template_string(INDEX_HTML)


@app.route("/analyze", methods=["POST"])
def analyze():
    try:
        ref_file = request.files.get("ref_file")
        sample_file = request.files.get("sample_file")
        ref_text = request.form.get("ref_text", "")
        sample_text = request.form.get("sample_text", "")

        ref_seq = parse_fasta_from_upload(ref_file, ref_text)
        sample_seq = parse_fasta_from_upload(sample_file, sample_text)
        if ref_seq is None or sample_seq is None:
            return jsonify({"error": "Reference or sample sequence missing."})

        variants, subs, indels = align_and_compare(ref_seq, sample_seq)
        ref_len = len(str(ref_seq))
        score = compute_risk_score(variants, subs, indels, ref_len)

        summary = {
            "total_variants": len(variants),
            "substitutions": subs,
            "indels": indels,
            "reference_length": ref_len,
            "risk_score": score,
        }

        return jsonify({"variants": variants, "summary": summary})
    except Exception as e:
        return jsonify({"error": str(e)})


if __name__ == "__main__":
    app.run(debug=True)