/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics.Utilities;

namespace Trinity
{
    public class Precursors : GraphML_List<Precursor>
    {
        public Precursors() { }
        private FDRizer<Precursor> uptimizer;

        public void OptimizePSMScoreRatios(DBOptions options, double desired_fdr)
        {
            List<Precursor> sortedPrecursorPrecision = new List<Precursor>(this);
            sortedPrecursorPrecision.Sort(Precursor.ComparePrecursorScore);
            double ratioPrecursorPrecision = FDRizer<Precursor>.ComputeAtFDR(sortedPrecursorPrecision, desired_fdr).Count / (double)this.Count;
           
            List<Precursor> sortedFragments = new List<Precursor>(this);
            sortedFragments.Sort(Precursor.CompareMatchingProductsFraction);
            double ratioFragments = FDRizer<Precursor>.ComputeAtFDR(sortedFragments, desired_fdr).Count / (double)this.Count;

            List<Precursor> sortedIntensities = new List<Precursor>(this);
            sortedIntensities.Sort(Precursor.CompareMatchingIntensityFraction);
            double ratioIntensities = FDRizer<Precursor>.ComputeAtFDR(sortedIntensities, desired_fdr).Count / (double)this.Count;

            List<Precursor> sortedPeptides = new List<Precursor>(this);
            sortedPeptides.Sort(Precursor.ComparePeptideScore);
            double ratioPeptideScore = FDRizer<Precursor>.ComputeAtFDR(sortedPeptides, desired_fdr).Count / (double)this.Count;

            double totalRatios = ratioPrecursorPrecision + ratioFragments + ratioIntensities + ratioPeptideScore;
            options.dPrecursor                  = ratioPrecursorPrecision   / totalRatios;
            options.dMatchingProductFraction    = ratioFragments            / totalRatios;
            options.dIntensityFraction          = ratioIntensities          / totalRatios;
            options.dPeptideScore               = ratioPeptideScore         / totalRatios;
            Console.WriteLine("New score ratios   ------------------------------------- ");
            Console.WriteLine("    PeptideSpectrumMatch.dPrecursor:                     " + options.dPrecursor);
            Console.WriteLine("    PeptideSpectrumMatch.dMatchingProductFraction:       " + options.dMatchingProductFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dIntensityFraction:             " + options.dIntensityFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dPeptideScore:                  " + options.dPeptideScore);
            Console.WriteLine("-------------------------------------------------------- ");
        }

        public List<Precursor> ComputeAtFDR(double desired_fdr, bool displayValues = true)
        {
            if (uptimizer == null)
            {
                List<Comparison<Precursor>> sorts = new List<Comparison<Precursor>>();
                sorts.Add(Precursor.OptimizedBestPSMComparison);
                //sorts.Add(Precursor.CompareOptimizedScore);
                sorts.Add(Precursor.CompareMatchingProducts);
                //sorts.Add(Precursor.CompareMatchingProductsFraction);
                sorts.Add(Precursor.ComparePrecursorScore);
                sorts.Add(Precursor.CompareMatchingIntensityFraction);
                sorts.Add(Precursor.ComparePeptideScore);
                sorts.Add(Precursor.CompareProteinScore);
                //sorts.Add(Precursor.CompareProductScore);
                //sorts.Add(Precursor.CompareScore);
                //sorts.Add(Precursor.CompareNbChargedPrecursor);//*/
                uptimizer = new FDRizer<Precursor>(this, sorts, null);
            }
            else
                uptimizer.ReStart();

            List<Precursor> fdrList = uptimizer.Launch(desired_fdr, displayValues);
            List<Precursor> sortedProbability = new List<Precursor>(this);
            sortedProbability.Sort(Precursor.CompareProbabilityScore);
            sortedProbability = FDRizer<Precursor>.ComputeAtFDR(sortedProbability, desired_fdr);
            if (sortedProbability.Count > fdrList.Count)
                fdrList = sortedProbability;

            return fdrList;
        }
    }

    public class Precursor : GraphML_Node, ITargetDecoy
    {
        public bool Decoy
        {
            get
            {
                PeptideSpectrumMatch psm = OptimizedBestPsm();
                if (psm != null)
                    return psm.Decoy;
                else
                    return true;
            }
        }

        public bool Target
        {
            get { return !Decoy; }
        }

        //MHC
        /*
        public static double dProduct = 0;
        public static double dPrecursor = 0;
        public static double dMatching = 0;//0.81;
        public static double dIntensity = 0;//0.61;
        public static double dProtein = 0;//0.00001;
        public static double dInvented = 0;//0;
        public static double dIsotope = 0;//0;
        public static double dScore = 0.5;//0;
        public static double dPeptideScore = 0.5;//0;//*/

        //Yeast
        /*
        public static double dProduct = 0.81;
        public static double dPrecursor = 0.21;
        public static double dMatching = 0.41;
        public static double dIntensity = 0.81;
        public static double dProtein = 0.0;
        public static double dInvented = 0;
        public static double dIsotope = 0;
        public static double dScore = 0;
        public static double dPeptideScore = 0;//*/

        public static int CompareMatchingIntensityFraction(Precursor left, Precursor right)
        {
            return -left.BestPSMMatchingIntensityFraction().CompareTo(right.BestPSMMatchingIntensityFraction());
        }
        public static int CompareProductScore(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().ProductScore.CompareTo(right.OptimizedBestPsm().ProductScore);
        }
        public static int ComparePrecursorScore(Precursor left, Precursor right)
        {
            return -left.BestPSMPrecursorScore().CompareTo(right.BestPSMPrecursorScore());
        }
        public static int CompareMatchingProductsFraction(Precursor left, Precursor right)
        {
            return -left.BestPSMMatchingProductsFraction().CompareTo(right.BestPSMMatchingProductsFraction());
        }
        public static int CompareMatchingProducts(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().MatchingProducts.CompareTo(right.OptimizedBestPsm().MatchingProducts);
        }
        public static int CompareProteinScore(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().ProteinScore.CompareTo(right.OptimizedBestPsm().ProteinScore);
        }
        public static int CompareProbabilityScore(Precursor left, Precursor right)
        {
            return -left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }
        public static int ComparePeptideScore(Precursor left, Precursor right)
        {
            return -left.BestPSMPeptideScore().CompareTo(right.BestPSMPeptideScore());
        }
        public static int CompareNbChargedPrecursor(Precursor left, Precursor right)
        {
            return -left.OtherCharges.Count.CompareTo(right.OtherCharges.Count);
        }
        public static int OptimizedBestPSMComparison(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsmScore().CompareTo(right.OptimizedBestPsmScore());
        }

        public double ProbabilityScore(Peptide peptide = null, bool checkMods = false)
        {/*
            double score = 0;
            foreach(PeptideSpectrumMatch psm in psms)
                if(
            foreach (PeptideSpectrumMatch psm in OptimizedBestPsms(peptide))
                score += (1 - score) * psm.OptimizedScore();
            return score;
            //*/
            return OptimizedBestPsmScore(peptide, checkMods);//TODO Reactivate Optimized Precursor Score
            //PeptideSpectrumMatch psmLeft = OptimizedBestPsm(peptide, checkMods);
            //if (psmLeft != null)
            //    return psmLeft.OptimizedScore();
            /*         (dIntensity * psmLeft.MatchingIntensityFraction +
                                dProduct * psmLeft.ProductScore +
                                dPrecursor * psmLeft.PrecursorScore +
                                dMatching * psmLeft.MatchingProductsFraction +
                                dProtein * psmLeft.ProteinScore +
                                dScore * psmLeft.MaxQuantScore() +
                                (Isotopes.Count > 0 ? dIsotope : 0) +
                                (Track.Invented ? dInvented : 0) +
                                dPeptideScore * psmLeft.PeptideScore);//*/
            //else
            //    return 0;//*/
        }

        public IEnumerable<PeptideSpectrumMatch> OptimizedBestPsms()
        {
            if (psms.Count > 0)
            {
                PeptideSpectrumMatch bestpsm = OptimizedBestPsm();
                double score = bestpsm.ProbabilityScore();
                for (int i = 0; i < psms.Count; i++)
                    if (psms[i].ProbabilityScore() >= score)
                        yield return psms[i];
            }
        }

        public static int OptimizedScoreComparison(Precursor left, Precursor right)
        {
            return -left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }

        public static int COMPTEUR = 0;
        //public GraphML_List<PeptideSpectrumMatch> psms;
        public PeptideSpectrumMatches psms;
        public PeptideSpectrumMatches psms_AllPossibilities;
        public Track Track;
        //public double Mz;
        public double Mass;
        //public double Rt;
        //public double Intensity;
        public int Charge;
        public GraphML_List<Precursor> Isotopes;
        public GraphML_List<Precursor> OtherCharges;
        public double MassShift;
        public int INDEX;
        public Sample sample;
        public Precursor()
        {
            this.psms_AllPossibilities = new PeptideSpectrumMatches();// GraphML_List<PeptideSpectrumMatch>();
            this.psms = new PeptideSpectrumMatches();// GraphML_List<PeptideSpectrumMatch>();
            Isotopes = new GraphML_List<Precursor>();
            OtherCharges = new GraphML_List<Precursor>();
        }
        public Precursor(Track track, int charge, Sample entry, double massShift = 0, GraphML_List<Precursor> isotopes = null)//double mz, double intensity, int charge = -1, double rt = 0, double massShift = 0, List<Precursor> isotopes = null)
        {
            INDEX = COMPTEUR++;
            this.Track = track;
            this.Charge = charge;
            this.sample = entry;

            this.psms_AllPossibilities = new PeptideSpectrumMatches();// GraphML_List<PeptideSpectrumMatch>();
            this.psms                  = new PeptideSpectrumMatches();// GraphML_List<PeptideSpectrumMatch>();
            //  this.Mz         = mz;
            //  this.Intensity  = intensity;
            //  this.Charge     = charge;
            this.Mass = Numerics.MassFromMZ(track.MZ, charge);
            //  this.Rt         = rt;           
            MassShift = massShift;
            if (isotopes == null)
                Isotopes = new GraphML_List<Precursor>();
            else
                Isotopes = isotopes;
            OtherCharges = new GraphML_List<Precursor>();
        }

        public int GetMostIntenseCharge()
        {
            double intensity = 0;
            int charge = -1;
            foreach (Precursor precursor in OtherCharges)
                if (precursor.Track.INTENSITY > intensity)
                {
                    intensity = precursor.Track.INTENSITY;
                    charge = precursor.Charge;
                }
            return charge;
        }
        

        public double OptimizedBestPsmScore(Peptide peptide = null, bool checkMods = false)
        {
            PeptideSpectrumMatch psm = OptimizedBestPsm(peptide, checkMods);
            if (psm != null)
                return psm.ProbabilityScore();
            else
                return 0;
        }

        public PeptideSpectrumMatch OptimizedBestPsm(Peptide peptide = null, bool checkMods = false)
        {
            if (psms.Count > 0)
            {
                if (peptide == null)
                {
                    double bestScore = psms[0].ProbabilityScore();
                    int bestProteinIndex = 0;
                    double tmpScore;
                    for (int i = 1; i < psms.Count; i++)
                    {
                        tmpScore = psms[i].ProbabilityScore();
                        if (tmpScore > bestScore || (tmpScore == bestScore && psms[i].Target))//.ProteinScore >= proteinScore))
                        {
                            bestProteinIndex = i;
                            bestScore = tmpScore;
                        }
                    }
                    return psms[bestProteinIndex];
                }
                else
                {
                    double bestScore = 0;
                    int bestProteinIndex = 0;
                    double proteinScore = 0;
                    for (int i = 0; i < psms.Count; i++)
                        if (peptide.IsSamePeptide(psms[i].Peptide, checkMods) && (psms[i].ProbabilityScore() > bestScore || (psms[i].ProbabilityScore() == bestScore && psms[i].Target)))//.ProteinScore > proteinScore)))
                        {
                            bestProteinIndex = i;
                            bestScore = psms[i].ProbabilityScore();
                            proteinScore = psms[i].ProteinScore;
                        }
                    if (bestScore > 0)
                        return psms[bestProteinIndex];
                }
            }
            return null;
        }

        public double BestPSMMatchingProductsFraction()
        {
            if (psms.Count > 0)
            {
                double bestScore = psms[0].MatchingProductsFraction;
                double tmpScore;
                for (int i = 1; i < psms.Count; i++)
                {
                    tmpScore = psms[i].MatchingProductsFraction;
                    if (tmpScore > bestScore)
                        bestScore = tmpScore;
                }
                return bestScore;
            }
            return 0.0;
        }
        public double BestPSMMatchingIntensityFraction()
        {
            if (psms.Count > 0)
            {
                double bestScore = psms[0].MatchingIntensityFraction;
                double tmpScore;
                for (int i = 1; i < psms.Count; i++)
                {
                    tmpScore = psms[i].MatchingIntensityFraction;
                    if (tmpScore > bestScore)
                        bestScore = tmpScore;
                }
                return bestScore;
            }
            return 0.0;
        }
        public double BestPSMPeptideScore()
        {
            if (psms.Count > 0)
            {
                double bestScore = psms[0].PeptideScore;
                double tmpScore;
                for (int i = 1; i < psms.Count; i++)
                {
                    tmpScore = psms[i].PeptideScore;
                    if (tmpScore > bestScore)
                        bestScore = tmpScore;
                }
                return bestScore;
            }
            return 0.0;
        }
        public double BestPSMPrecursorScore()
        {
            if (psms.Count > 0)
            {
                double bestScore = psms[0].PrecursorScore;
                double tmpScore;
                for (int i = 1; i < psms.Count; i++)
                {
                    tmpScore = psms[i].PrecursorScore;
                    if (tmpScore > bestScore)
                        bestScore = tmpScore;
                }
                return bestScore;
            }
            return 0.0;
        }

    }
}
