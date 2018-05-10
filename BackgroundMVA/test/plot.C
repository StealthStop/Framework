#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>

//This is a helper function which will keep the plot from overlapping with the legend
void smartMax(const TH1 * const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax, const bool error)
{
    const bool isLog = p->GetLogy();
    double min = 9e99;
    double max = -9e99;
    double pThreshMax = -9e99;
    int threshold = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        double bin = 0.0;
        if(error) bin = h->GetBinContent(i) + h->GetBinError(i);
        else      bin = h->GetBinContent(i);
        if(bin > max) max = bin;
        else if(bin > 1e-10 && bin < min) min = bin;
        if(i >= threshold && bin > pThreshMax) pThreshMax = bin;
    }

    gpThreshMax = std::max(gpThreshMax, pThreshMax);
    gmax = std::max(gmax, max);
    gmin = std::min(gmin, min);
}

//Class to hold TH1* with various helper functions 
class histInfo
{
public:
    std::string legEntry, histFile, histName, drawOptions;
    int color, rebin;
    std::shared_ptr<TH1> h;

    //helper function to get histogram from file and configure its optional settings
    void retrieveHistogram()
    {
        //Open the file for this histogram
        TFile *f = TFile::Open(histFile.c_str());

        //check that the file was opened successfully
        if(!f)
        {
            printf("File \"%s\" could not be opened!!!\n", histFile.c_str());
            h = nullptr;
            return;
        }

        //get the histogram from the file
        h.reset(static_cast<TH1*>(f->Get(histName.c_str())));

        //with the histogram retrieved, close the file
        f->Close();
        delete f;

        //check that the histogram was retireved from the file successfully
        if(!h)
        {
            printf("Histogram \"%s\" could not be found in file \"%s\"!!!\n", histName.c_str(), histFile.c_str());
            return;
        }

        //set the histogram color
        h->SetLineColor(color);
        h->SetLineWidth(3);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(20);

        // rebin the histogram if desired
        if(rebin > 0) h->Rebin(rebin);
    }

    //helper function for axes
    void setupAxes(double xOffset, double yOffset, double xTitle, double yTitle, double xLabel, double yLabel)
    {
        h->SetStats(0);
        h->SetTitle(0);
        h->GetXaxis()->SetTitleOffset(xOffset);
        h->GetYaxis()->SetTitleOffset(yOffset);
        h->GetXaxis()->SetTitleSize(xTitle);
        h->GetYaxis()->SetTitleSize(yTitle);
        h->GetXaxis()->SetLabelSize(xLabel);
        h->GetYaxis()->SetLabelSize(yLabel);
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) h->GetXaxis()->SetNdivisions(6, 5, 0);
    }

    //helper function for pads
    void setupPad(double left, double right, double top, double bottom)
    {
        gPad->SetLeftMargin(left);
        gPad->SetRightMargin(right);
        gPad->SetTopMargin(top);
        gPad->SetBottomMargin(bottom);
        gPad->SetTicks(1,1);
    }

    //wrapper to draw histogram
    void draw(const std::string& additionalOptions = "", bool noSame = false) const
    {
        h->Draw(((noSame?"":"same " + drawOptions + " " + additionalOptions)).c_str());
    }

    void setFillColor(int newColor = -1)
    {
        if(newColor >= 0) h->SetFillColor(newColor);
        else              h->SetFillColor(color);
    }

    void normalize(double num = 1.0)
    {
        h->Scale(num/h->Integral());
    }

    histInfo(const std::string& legEntry, const std::string& histFile, const std::string& drawOptions, const int color) : legEntry(legEntry), histFile(histFile), histName(""), drawOptions(drawOptions), color(color), rebin(-1), h(nullptr)
    {
    }

    histInfo(const std::string& legEntry, const std::string& histFile, const std::string& histName, const std::string& drawOptions, const int color) : legEntry(legEntry), histFile(histFile), histName(histName), drawOptions(drawOptions), color(color), rebin(-1), h(nullptr)
    {
    }

    histInfo(TH1* h) : legEntry(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(0), rebin(0), h(h)
    {
    }

    ~histInfo()
    {
    }
};

class Plotter
{
private:
    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> fisherG_;
    
public:
    Plotter(std::vector<histInfo>&  fisherG) : fisherG_(fisherG) {}
    Plotter(std::vector<histInfo>&& fisherG) : fisherG_(fisherG) {}

    void plot(const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        //switch to the canvas to ensure it is the active object
        c->cd();

        //Create TLegend: TLegend(x1, y1, x2, y2)
        TLegend *leg = new TLegend(0.50, 0.76, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        //signal 
        for(auto& entry : fisherG_)
        {
            //get new histogram
            entry.rebin = rebin;
            entry.retrieveHistogram();
            entry.normalize();

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L");
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        }

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, fisherG_[0].h->GetBinLowEdge(1), fisherG_[0].h->GetBinLowEdge(fisherG_[0].h->GetNbinsX()) + fisherG_[0].h->GetBinWidth(fisherG_[0].h->GetNbinsX())));
        // set pad margins: setupPad(left, right, top, bottom)
        dummy.setupPad(0.17, 0.06, 0.08, 0.17);
        dummy.setupAxes(1.0, 1.4, 0.06, 0.06, 0.05, 0.05);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTickLength(0.03);
        dummy.h->GetYaxis()->SetTickLength(0.03);

        //Set the y-range of the histogram
        if(isLogY)
        {
            double locMin = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin = legSpan + log10(locMin);
            if(log10(lmax) > legMin)
            {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
            }
            dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
        }
        else
        {
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
            dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.3);
        }
        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //plot signal histograms
        for(const auto& entry : fisherG_)
        {
            entry.draw();
        }

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        c->Print( ("outputPlots/"+xAxisLabel+".png").c_str() );

        //clean up dynamic memory
        delete c;
        delete leg;
    }
};

int main()
{
    const std::string& file = "outputfiles/tmva-train-example-output.root";
    const std::string& path = "fisherLoader/InputVariables_Id/";

    const std::vector<std::string>& plotTypes = {"fwm2_top6", "fwm3_top6", "fwm4_top6", "fwm5_top6", "fwm6_top6", "jmt_ev0_top6", "jmt_ev1_top6", "jmt_ev2_top6"};
        
    for(const auto& t : plotTypes)
    {
        std::vector<histInfo> fisherG = {
            {"Signal",     file, path+t+"__Signal_Id",     "hist", kBlack},
            {"Background", file, path+t+"__Background_Id", "hist", kRed},
        };
        Plotter plt(fisherG);
        plt.plot(t, "A.U.");
    }
}
