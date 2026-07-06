/* TCR Germline Gene Mapper — frontend logic (vanilla JS, no build) */

(() => {
  "use strict";

  // Detect API base URL: when served by Vercel, same origin. When running
  // locally with `python -m http.server` in public/ on port 5173, the API
  // lives on port 8000. Allow override via ?api= query param.
  const params = new URLSearchParams(location.search);
  const API_BASE = params.get("api") ||
    (location.protocol === "http:" && location.port === "5173"
      ? "http://localhost:8000"
      : "");

  const $ = (id) => document.getElementById(id);

  const state = {
    file: null,        // File object
    fileBytes: null,   // ArrayBuffer
    uploadGuess: null, // server's guess
    lastResult: null,  // last process result
  };

  // ---------- DOM refs ----------
  const els = {
    dropzone: $("dropzone"),
    fileInput: $("file-input"),
    uploadStatus: $("upload-status"),
    statusFilename: $("status-filename"),
    statusSize: $("status-size"),
    statusFormat: $("status-format"),
    statusMolecule: $("status-molecule"),
    statusNotes: $("status-notes"),

    stepUpload: $("step-upload"),
    stepConfirm: $("step-confirm"),
    stepResults: $("step-results"),

    confirmFormat: $("confirm-format"),
    confirmMolecule: $("confirm-molecule"),
    confirmChains: $("confirm-chains"),
    moleculeRow: $("molecule-row"),
    chainsRow: $("chains-row"),

    btnBack: $("btn-back"),
    btnProcess: $("btn-process"),
    btnReset: $("btn-reset"),
    btnDownloadJson: $("btn-download-json"),

    processingIndicator: $("processing-indicator"),
    resultsArea: $("results-area"),
    resultsElapsed: $("results-elapsed"),
    resultsChaincount: $("results-chaincount"),
    resultsChains: $("results-chains"),
    warningsDetails: $("warnings-details"),
    warningsList: $("warnings-list"),
    warningCount: $("warning-count"),
    referenceDetails: $("reference-details"),
    referencePre: $("reference-pre"),

    apiDocsLink: $("api-docs-link"),
    apiDialog: $("api-dialog"),
    apiDialogClose: $("api-dialog-close"),
  };

  // ---------- Utility ----------
  function formatBytes(n) {
    if (n < 1024) return `${n} B`;
    if (n < 1024 * 1024) return `${(n / 1024).toFixed(1)} KB`;
    return `${(n / (1024 * 1024)).toFixed(2)} MB`;
  }

  function showStep(step) {
    [els.stepUpload, els.stepConfirm, els.stepResults].forEach((el) => {
      el.hidden = true;
    });
    step.hidden = false;
    step.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function setError(msg) {
    // Lightweight inline error display
    let errEl = document.querySelector(".global-error");
    if (!errEl) {
      errEl = document.createElement("div");
      errEl.className = "card global-error";
      errEl.style.borderColor = "var(--bad)";
      errEl.style.color = "#fecaca";
      document.querySelector("main").prepend(errEl);
    }
    errEl.textContent = msg;
    errEl.scrollIntoView({ behavior: "smooth", block: "start" });
    setTimeout(() => errEl.remove(), 8000);
  }

  // ---------- Step 1: file upload ----------
  ["dragenter", "dragover"].forEach((ev) => {
    els.dropzone.addEventListener(ev, (e) => {
      e.preventDefault();
      e.stopPropagation();
      els.dropzone.classList.add("drag-over");
    });
  });
  ["dragleave", "drop"].forEach((ev) => {
    els.dropzone.addEventListener(ev, (e) => {
      e.preventDefault();
      e.stopPropagation();
      els.dropzone.classList.remove("drag-over");
    });
  });
  els.dropzone.addEventListener("drop", (e) => {
    const files = e.dataTransfer.files;
    if (files && files.length) handleFileSelected(files[0]);
  });
  els.dropzone.addEventListener("click", () => els.fileInput.click());
  els.dropzone.addEventListener("keydown", (e) => {
    if (e.key === "Enter" || e.key === " ") {
      e.preventDefault();
      els.fileInput.click();
    }
  });
  els.fileInput.addEventListener("change", (e) => {
    if (e.target.files && e.target.files.length) {
      handleFileSelected(e.target.files[0]);
    }
  });

  async function handleFileSelected(file) {
    state.file = file;

    // Reset prior results
    state.lastResult = null;
    els.resultsArea.hidden = true;
    els.processingIndicator.hidden = true;

    if (file.size > 4_500_000) {
      setError(`File too large: ${formatBytes(file.size)}. Max is 4.5 MB (Vercel Hobby tier).`);
      return;
    }

    // Show filename immediately
    els.statusFilename.textContent = file.name;
    els.statusSize.textContent = formatBytes(file.size);
    els.uploadStatus.hidden = false;

    // Send to /api/upload for sniffing
    try {
      const fd = new FormData();
      fd.append("file", file);
      const resp = await fetch(`${API_BASE}/api/upload`, { method: "POST", body: fd });
      if (!resp.ok) {
        const err = await resp.json().catch(() => ({ error: resp.statusText }));
        throw new Error(err.error || `Upload failed (${resp.status})`);
      }
      const data = await resp.json();
      state.uploadGuess = data;

      els.statusFormat.textContent = data.file_format || "unknown";
      els.statusMolecule.textContent = data.suggested_molecule || "unknown";

      els.statusNotes.innerHTML = "";
      (data.sniffed_notes || []).forEach((n) => {
        const li = document.createElement("li");
        li.textContent = n;
        els.statusNotes.appendChild(li);
      });

      // Pre-fill step 2 with the guess
      if (data.file_format === "pdb" || data.file_format === "mmcif") {
        els.confirmFormat.value = data.file_format;
        els.confirmMolecule.value =
          data.suggested_molecule === "pmhc_tcr_complex"
            ? "pmhc_tcr_complex"
            : "tcr_only";
      } else if (data.file_format === "fasta") {
        els.confirmFormat.value = "fasta";
        els.confirmMolecule.value = "fasta";
      }
      updateMoleculeRowVisibility();

      // Move to step 2
      els.stepConfirm.hidden = false;
      els.stepConfirm.scrollIntoView({ behavior: "smooth", block: "start" });
    } catch (err) {
      setError(`Upload failed: ${err.message}`);
    }
  }

  // ---------- Step 2: confirm + dispatch to step 3 ----------
  els.confirmFormat.addEventListener("change", updateMoleculeRowVisibility);
  els.confirmMolecule.addEventListener("change", updateMoleculeRowVisibility);

  function updateMoleculeRowVisibility() {
    const fmt = els.confirmFormat.value;
    const mol = els.confirmMolecule.value;
    // Molecule selector only meaningful for PDB/mmCIF
    els.moleculeRow.hidden = !(fmt === "pdb" || fmt === "mmcif");
    // Chain IDs only when pMHC complex (or auto with PDB)
    els.chainsRow.hidden = !(fmt === "pdb" || fmt === "mmcif") ||
      !(mol === "pmhc_tcr_complex" || mol === "auto");
  }

  els.btnBack.addEventListener("click", () => {
    showStep(els.stepUpload);
  });

  els.btnProcess.addEventListener("click", async () => {
    if (!state.file) {
      setError("No file selected.");
      return;
    }
    showStep(els.stepResults);
    els.processingIndicator.hidden = false;
    els.resultsArea.hidden = true;

    const fd = new FormData();
    fd.append("file", state.file);
    fd.append("file_format", els.confirmFormat.value);
    fd.append("molecule", els.confirmMolecule.value);
    if (els.confirmChains.value.trim()) {
      fd.append("tcr_chain_ids", els.confirmChains.value.trim());
    }

    try {
      const resp = await fetch(`${API_BASE}/api/process`, { method: "POST", body: fd });
      if (!resp.ok) {
        const err = await resp.json().catch(() => ({ error: resp.statusText }));
        throw new Error(err.error || `Processing failed (${resp.status})`);
      }
      const data = await resp.json();
      state.lastResult = data;
      renderResults(data);
    } catch (err) {
      setError(`Processing failed: ${err.message}`);
      els.processingIndicator.hidden = true;
      els.resultsArea.hidden = false;
      els.resultsChains.innerHTML = `<p style="color:#fecaca">Pipeline error: ${err.message}</p>`;
    }
  });

  // ---------- Step 3: render results ----------
  function renderResults(data) {
    els.processingIndicator.hidden = true;
    els.resultsArea.hidden = false;

    els.resultsElapsed.textContent = `Elapsed: ${data.elapsed_ms} ms`;
    els.resultsChaincount.textContent = `Chains: ${data.chains.length}`;

    // Chains
    els.resultsChains.innerHTML = "";
    if (!data.chains.length) {
      els.resultsChains.innerHTML = `<p class="muted">No TCR chains found.</p>`;
    } else {
      data.chains.forEach((chain) => {
        els.resultsChains.appendChild(renderChain(chain));
      });
    }

    // Warnings
    const warnings = data.warnings || [];
    els.warningCount.textContent = warnings.length;
    els.warningsList.innerHTML = "";
    warnings.forEach((w) => {
      const li = document.createElement("li");
      li.textContent = w;
      els.warningsList.appendChild(li);
    });
    els.warningsDetails.hidden = warnings.length === 0;

    // Reference info
    els.referencePre.textContent = JSON.stringify(data.reference_info, null, 2);
  }

  function renderChain(chain) {
    const div = document.createElement("div");
    div.className = "chain-result";

    const ctype = (chain.chain_type || "unknown").toLowerCase();
    const geneCells = [
      { label: "V gene", value: chain.v_allele, gene: chain.v_gene, ident: chain.v_identity_pct },
      { label: "J gene", value: chain.j_allele, gene: chain.j_gene, ident: chain.j_identity_pct },
      { label: "C gene", value: chain.c_allele, gene: chain.c_gene, ident: chain.c_identity_pct },
    ];
    if (chain.d_allele) {
      geneCells.push({
        label: "D gene",
        value: chain.d_allele,
        gene: chain.d_gene,
        ident: chain.d_identity_pct,
        lowConf: chain.d_confidence,
      });
    }

    div.innerHTML = `
      <h4>
        <span>Chain ${escapeHtml(chain.chain_id)}</span>
        <span class="chain-type-badge ${ctype}">${escapeHtml(chain.chain_type || "unknown")}</span>
      </h4>
      <div class="gene-grid">
        ${geneCells.map(renderGeneCell).join("")}
      </div>
      <div class="cdr3-block">
        <div class="label">Inferred CDR3</div>
        <div class="seq ${chain.cdr3_sequence ? "" : "empty"}">
          ${chain.cdr3_sequence ? escapeHtml(chain.cdr3_sequence) : "(not extracted)"}
        </div>
      </div>
      ${(chain.notes && chain.notes.length) ? `
        <ul class="chain-notes">
          ${chain.notes.map((n) => `<li>${escapeHtml(n)}</li>`).join("")}
        </ul>
      ` : ""}
    `;
    return div;
  }

  function renderGeneCell(cell) {
    const val = cell.value || "—";
    const genePart = cell.gene && cell.gene !== cell.value ? ` <span class="muted-val">(${escapeHtml(cell.gene)})</span>` : "";
    const ident = typeof cell.ident === "number" && cell.ident > 0 ? cell.ident.toFixed(1) + "%" : "—";
    const barPct = typeof cell.ident === "number" ? Math.max(0, Math.min(100, cell.ident)) : 0;
    const lowConf = cell.lowConf ? `<span class="d-conf-low">${escapeHtml(cell.lowConf)}</span>` : "";
    return `
      <div class="gene-cell">
        <div class="label">${escapeHtml(cell.label)}</div>
        <div class="value">${escapeHtml(val)}${genePart}</div>
        <span class="ident">identity: ${ident}</span>
        <span class="ident-bar"><span style="width:${barPct}%"></span></span>
        ${lowConf}
      </div>
    `;
  }

  function escapeHtml(s) {
    if (s === null || s === undefined) return "";
    return String(s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  // ---------- Reset / download ----------
  els.btnReset.addEventListener("click", () => {
    state.file = null;
    state.uploadGuess = null;
    state.lastResult = null;
    els.fileInput.value = "";
    els.uploadStatus.hidden = true;
    els.resultsArea.hidden = true;
    els.processingIndicator.hidden = true;
    showStep(els.stepUpload);
  });

  els.btnDownloadJson.addEventListener("click", () => {
    if (!state.lastResult) return;
    const blob = new Blob([JSON.stringify(state.lastResult, null, 2)], {
      type: "application/json",
    });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "tcr_germline_mapping.json";
    a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });

  // ---------- API docs dialog ----------
  els.apiDocsLink.addEventListener("click", (e) => {
    e.preventDefault();
    els.apiDialog.showModal();
  });
  els.apiDialogClose.addEventListener("click", () => els.apiDialog.close());
  els.apiDialog.addEventListener("click", (e) => {
    if (e.target === els.apiDialog) els.apiDialog.close();
  });

  // ---------- Health check on load (non-blocking) ----------
  fetch(`${API_BASE}/api/health`)
    .then((r) => r.json())
    .then((data) => {
      console.log("[health]", data);
    })
    .catch((err) => {
      console.warn("[health] API unreachable:", err.message);
    });
})();
