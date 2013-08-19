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

        public void OptimizePSMScoreRatios(double desired_fdr)
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
            PeptideSpectrumMatch.dPrecursor                 = ratioPrecursorPrecision   / totalRatios;
            PeptideSpectrumMatch.dMatchingProductFraction   = ratioFragments            / totalRatios;
            PeptideSpectrumMatch.dIntensityFraction         = ratioIntensities          / totalRatios;
            PeptideSpectrumMatch.dPeptideScore              = ratioPeptideScore         / totalRatios;
            Console.WriteLine("New score ratios   ------------------------------------- ");
            Console.WriteLine("    PeptideSpectrumMatch.dPrecursor: " + PeptideSpectrumMatch.dPrecursor);
            Console.WriteLine("    PeptideSpectrumMatch.dMatchingProductFraction: " + PeptideSpectrumMatch.dMatchingProductFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dIntensityFraction: " + PeptideSpectrumMatch.dIntensityFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dPeptideScore: " + PeptideSpectrumMatch.dPeptideScore);
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
            sortedProbability.Sort(Precursor.CompareOptimizedScore);
            sortedProbability = FDRizer<Precursor>.ComputeAtFDR(sortedProbability, desired_fdr);
            if (sortedProbability.Count > fdrList.Count)
                fdrList = sortedProbability;

            return fdrList;
        }
    }

    public class Precursor : GraphML_Node, ITargetDecoy
    {
        public double Score
        {
            get { return ScoreFct(); }
        }

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
            return -left.OptimizedBestPsm().MatchingIntensityFraction.CompareTo(right.OptimizedBestPsm().MatchingIntensityFraction);
        }
        public static int CompareProductScore(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().ProductScore.CompareTo(right.OptimizedBestPsm().ProductScore);
        }
        public static int ComparePrecursorScore(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().PrecursorScore.CompareTo(right.OptimizedBestPsm().PrecursorScore);
        }
        public static int CompareMatchingProductsFraction(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().MatchingProductsFraction.CompareTo(right.OptimizedBestPsm().MatchingProductsFraction);
        }
        public static int CompareMatchingProducts(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().MatchingProducts.CompareTo(right.OptimizedBestPsm().MatchingProducts);
        }
        public static int CompareProteinScore(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().ProteinScore.CompareTo(right.OptimizedBestPsm().ProteinScore);
        }
        public static int CompareScore(Precursor left, Precursor right)
        {
            return -left.ScoreFct().CompareTo(right.ScoreFct());
        }
        public static int CompareOptimizedScore(Precursor left, Precursor right)
        {
            return -left.OptimizedScore().CompareTo(right.OptimizedScore());
        }
        public static int ComparePeptideScore(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsm().PeptideScore.CompareTo(right.OptimizedBestPsm().PeptideScore);
        }
        public static int CompareNbChargedPrecursor(Precursor left, Precursor right)
        {
            return -left.OtherCharges.Count.CompareTo(right.OtherCharges.Count);
        }
        public static int OptimizedBestPSMComparison(Precursor left, Precursor right)
        {
            return -left.OptimizedBestPsmScore().CompareTo(right.OptimizedBestPsmScore());
        }

        public double OptimizedScore(Peptide peptide = null, bool checkMods = false)
        {
            //return OptimizedBestPsmScore(peptide, false);//TODO Reactivate Optimized Precursor Score
            PeptideSpectrumMatch psmLeft = OptimizedBestPsm(peptide, checkMods);
            if (psmLeft != null)
                return psmLeft.OptimizedScore();/*         (dIntensity * psmLeft.MatchingIntensityFraction +
                                dProduct * psmLeft.ProductScore +
                                dPrecursor * psmLeft.PrecursorScore +
                                dMatching * psmLeft.MatchingProductsFraction +
                                dProtein * psmLeft.ProteinScore +
                                dScore * psmLeft.MaxQuantScore() +
                                (Isotopes.Count > 0 ? dIsotope : 0) +
                                (Track.Invented ? dInvented : 0) +
                                dPeptideScore * psmLeft.PeptideScore);//*/
            else
                return 0;//*/
        }
        public static int OptimizedScoreComparison(Precursor left, Precursor right)
        {
            return -left.OptimizedScore().CompareTo(right.OptimizedScore());
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

        public static int DescendingScoreComparison(Precursor left, Precursor right)
        {
            return -(left.ScoreFct(null).CompareTo(right.ScoreFct(null)));
        }

        public double ScoreFct(Peptide peptide = null, bool computeIsotopes = true)
        {
            return OptimizedScore(peptide);
        }
        
        public IEnumerable<PeptideSpectrumMatch> OptimizedBestPsms(Peptide peptide = null)
        {
            if (psms.Count > 0)
            {
                PeptideSpectrumMatch bestpsm = OptimizedBestPsm(peptide);
                double score = bestpsm.OptimizedScore();
                for (int i = 0; i < psms.Count; i++)
                    if (psms[i].OptimizedScore() >= score)
                        yield return psms[i];
            }
        }

        public double OptimizedBestPsmScore(Peptide peptide = null, bool checkMods = false)
        {
            PeptideSpectrumMatch psm = OptimizedBestPsm(peptide, checkMods);
            if (psm != null)
                return psm.OptimizedScore();
            else
                return 0;
        }

        public PeptideSpectrumMatch OptimizedBestPsm(Peptide peptide = null, bool checkMods = false)
        {
            if (psms.Count > 0)
            {
                if (peptide == null)
                {
                    double bestScore = psms[0].OptimizedScore();
                    int bestProteinIndex = 0;
                    double tmpScore;
                    for (int i = 1; i < psms.Count; i++)
                    {
                        tmpScore = psms[i].OptimizedScore();
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
                        if (peptide.IsSamePeptide(psms[i].Peptide, checkMods) && (psms[i].OptimizedScore() > bestScore  || (psms[i].OptimizedScore() == bestScore && psms[i].Target)))//.ProteinScore > proteinScore)))
                        {
                            bestProteinIndex = i;
                            bestScore = psms[i].OptimizedScore();
                            proteinScore = psms[i].ProteinScore;
                        }
                    if(bestScore > 0)
                        return psms[bestProteinIndex];
                }
            }
            return null;
        }
    }
}
