
// using  a weird default name so that code crash in cse a correct rootfile is not provided 
void PlotPulls(TString filename="pulls_none.root", TString outdir="", TString postfix_="",int numberOfCanvas=3){ 
  
  TString plotdir = outdir;
    TFile file(filename,"READ");
    TCanvas *c = (TCanvas*)file.Get("nuisances");
    c->ls(); //check inside the c canvas
    c->Size(1250,400);
    c->SetBottomMargin(0.35);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    //c->GetXaxis->SetRangeUser(0,5)
    
    TH1F *h1 = (TH1F*)c->GetPrimitive("prefit_nuisancs");
    h1->LabelsOption("v");
    //h1->SetAxisRange(0, 5, "X");
    //cout << "print " << endl;
    //cout << " testing " << h1->GetXaxis()->GetNbins()  << endl;
    int numberOfNuisance = h1->GetXaxis()->GetNbins();
    TLegend leg1 = TLegend(0.6, 0.74, 0.89, 0.89);
    TLegend *leg2 = (TLegend*)(c->GetPrimitive("TPave"));
    leg1.Copy(*leg2);
    
    TPaveText *pt = new TPaveText(0.0877181,0.9,0.9580537,0.96,"brNDC");
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetFillStyle(0);
    pt->SetTextFont(52);
    
    double cmstextSize = 0.07;
    double preliminarytextfize = cmstextSize * 0.7;
    double lumitextsize = cmstextSize *0.7;
    pt->SetTextSize(cmstextSize);
    pt->AddText(0.01,0.57,"#font[61]{CMS}");
    
    TPaveText *pt1 = new TPaveText(0.0877181,0.905,0.9580537,0.96,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetTextAlign(12);
    pt1->SetFillStyle(0);
    pt1->SetTextFont(52);
    pt1->SetTextSize(preliminarytextfize);
    //pt1->AddText(0.155,0.4,"Preliminary");
    pt1->AddText(0.125,0.4,"Internal");
    
    TPaveText *pt2 = new TPaveText(0.0877181,0.9,0.8280537,0.96,"brNDC");
    pt2->SetBorderSize(0);
    pt2->SetTextAlign(12);
    pt2->SetFillStyle(0);
    pt2->SetTextFont(42);
    pt2->SetTextSize(lumitextsize);
    pt2->AddText(0.81, 0.5, "41.5 fb^{-1} (13 TeV)");
    
    TPaveText *pt3 = new TPaveText(0.0377181,0.85,0.9580537,0.88,"brNDC");
    pt3->SetBorderSize(0);
    pt3->SetTextAlign(12);
    pt3->SetFillStyle(0);
    pt3->SetTextFont(42);
    pt3->SetTextSize(lumitextsize);
    TString latexText = "Data"; 
    
    if (filename.Contains("data")) latexText = "Data";
    if (filename.Contains("asimov")) latexText = "Asimov";
    
    pt3->AddText(0.1,0.4, latexText+" Fit");
    
    h1->Draw("same");
    leg2->Draw();
    pt->Draw();
    pt1->Draw();
    pt2->Draw();
    pt3->Draw();

    //int numberOfCanvas = 2 ; 
    int nuisanceRange = numberOfNuisance/numberOfCanvas;
    for (int i =1 ; i < numberOfCanvas+1 ; i++)
	{
	int dummy = 0;
	string postfix= to_string(i);
	string prefix= to_string(i-1);
	int last = nuisanceRange * i;
	int start = nuisanceRange * ( i -1 );
	cout << "start  "  << start << "  last  " << last << endl;
	if (i==1) {
		h1->SetAxisRange(start, last, "X");
		}
	if (i!=1 && i < numberOfCanvas) {
		h1->SetAxisRange(start+1, last, "X");
		}
	if (i == numberOfCanvas){
		h1->SetAxisRange(start+1, numberOfNuisance+1, "X");
		}
	c->Update();
        c->Modified();
	c->SaveAs(plotdir+filename.ReplaceAll(".root","_"+postfix+"_.pdf").ReplaceAll("_"+prefix+"_",""));
	c->SaveAs(plotdir+filename.ReplaceAll(".pdf",".png"));
	c->SaveAs(plotdir+filename.ReplaceAll(".png",".root"));
	
	}

    //h1->SetAxisRange(0, 15, "X");
   
    //c->Update();
    //c->Modified();
    
    //c->SaveAs(plotdir+filename.ReplaceAll(".root",".pdf"));
    //c->SaveAs(plotdir+filename.ReplaceAll(".pdf",".png"));
    //c->SaveAs(plotdir+filename.ReplaceAll(".png",".root"));
    //c->SaveAs(plotdir+filename.ReplaceAll(".root",".C"));


}
