void make_runbx_pdf() {
  // Open file
  TFile *f = TFile::Open("runbx_map.root");
  if (!f || f->IsZombie()) {
    Error("make_run_vs_bx", "Cannot open runbx_map.root");
    return;
  }

  // Get THnSparse
  THnSparse *h = (THnSparse*)f->Get("h2runbx");
  if (!h) {
    Error("make_run_vs_bx", "h2runbx not found");
    return;
  }

  // Optional: restrict run range (recommended!)
  // Example: one run
  // h->GetAxis(0)->SetRangeUser(395000, 395000);

  // Or a run window
  // h->GetAxis(0)->SetRangeUser(392000, 394000);

  // Project to 2D: X=run (axis 0), Y=BX (axis 1)
  TH2D *h2 = (TH2D*)h->Projection(1, 0);
  h2->SetTitle(";Run number;Bunch crossing");

  // Draw
  TCanvas *c = new TCanvas("c", "Run vs BX", 1200, 800);
  c->SetRightMargin(0.15);
  h2->Draw("COLZ");

  c->SaveAs("run_vs_bx.pdf");

  // Reset axis ranges (good practice)
  h->GetAxis(0)->SetRange(0, -1);
}