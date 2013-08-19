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
    public class PeptideSpectrumMatches : GraphML_List<PeptideSpectrumMatch>
    {
        FDRizer<PeptideSpectrumMatch> uptimizer;
        public PeptideSpectrumMatches() { }
        public PeptideSpectrumMatches(IEnumerable<PeptideSpectrumMatch> list) : base(list) { }

        public void OptimizePSMScoreRatios(double desired_fdr)
        {
            List<PeptideSpectrumMatch> sortedPrecursorPrecision = new List<PeptideSpectrumMatch>(this);
            sortedPrecursorPrecision.Sort(PeptideSpectrumMatches.ComparePrecursorScore);
            double ratioPrecursorPrecision = FDRizer<PeptideSpectrumMatch>.ComputeAtFDR(sortedPrecursorPrecision, desired_fdr).Count / (double)this.Count;

            List<PeptideSpectrumMatch> sortedFragments = new List<PeptideSpectrumMatch>(this);
            sortedFragments.Sort(PeptideSpectrumMatches.CompareMatchingProductsFraction);
            double ratioFragments = FDRizer<PeptideSpectrumMatch>.ComputeAtFDR(sortedFragments, desired_fdr).Count / (double)this.Count;

            List<PeptideSpectrumMatch> sortedIntensities = new List<PeptideSpectrumMatch>(this);
            sortedIntensities.Sort(PeptideSpectrumMatches.CompareMatchingIntensityFraction);
            double ratioIntensities = FDRizer<PeptideSpectrumMatch>.ComputeAtFDR(sortedIntensities, desired_fdr).Count / (double)this.Count;

            List<PeptideSpectrumMatch> sortedPeptides = new List<PeptideSpectrumMatch>(this);
            sortedPeptides.Sort(PeptideSpectrumMatches.ComparePeptideScore);
            double ratioPeptideScore = FDRizer<PeptideSpectrumMatch>.ComputeAtFDR(sortedPeptides, desired_fdr).Count / (double)this.Count;

            double totalRatios = ratioPrecursorPrecision + ratioFragments + ratioIntensities + ratioPeptideScore;
            PeptideSpectrumMatch.dPrecursor = ratioPrecursorPrecision / totalRatios;
            PeptideSpectrumMatch.dMatchingProductFraction = ratioFragments / totalRatios;
            PeptideSpectrumMatch.dIntensityFraction = ratioIntensities / totalRatios;
            PeptideSpectrumMatch.dPeptideScore = ratioPeptideScore / totalRatios;
            Console.WriteLine("New score ratios   ------------------------------------- ");
            Console.WriteLine("    PeptideSpectrumMatch.dPrecursor: " + PeptideSpectrumMatch.dPrecursor);
            Console.WriteLine("    PeptideSpectrumMatch.dMatchingProductFraction: " + PeptideSpectrumMatch.dMatchingProductFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dIntensityFraction: " + PeptideSpectrumMatch.dIntensityFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dPeptideScore: " + PeptideSpectrumMatch.dPeptideScore);
            Console.WriteLine("-------------------------------------------------------- ");
        }
        public static int CompareMatchingIntensityFraction(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingIntensityFraction.CompareTo(right.MatchingIntensityFraction);
        }
        public List<PeptideSpectrumMatch> ComputeAtFDR(double desired_fdr, bool precursorErrorOnly, bool display = false)
        {
            if (uptimizer == null)
            {
                List<Comparison<PeptideSpectrumMatch>> sorts = new List<Comparison<PeptideSpectrumMatch>>();
                if (precursorErrorOnly)
                    sorts.Add(CompareProductScore);
                else
                {
                    //sorts.Add(CompareMatchingIntensity);
                    sorts.Add(CompareMaxQuantScore);
                    //                sorts.Add(ComparePrecursorScore);
                    //                sorts.Add(CompareMatchingProductsFraction);
                    //                sorts.Add(CompareMatchingProducts);
                    //sorts.Add(CompareProteinScore);
                    //sorts.Add(ComparePeptideScore);
                }
                uptimizer = new FDRizer<PeptideSpectrumMatch>(this, sorts, null);
            }
            return uptimizer.Launch(desired_fdr, display);
        }
        /*
            return              dIntensity * MatchingIntensityFraction +
                                dProduct * ProductScore +
                                dPrecursor * PrecursorScore +
                                dMatchingProductFraction * MatchingProductsFraction +
                                dMatchingProduct * MatchingProducts +
                                dProtein * ProteinScore +
                                dPeptideScore * PeptideScore;*/
        public static int CompareTrapDistance(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.Query.TrapDistance.CompareTo(right.Query.TrapDistance);
        }
        public static int CompareMatchingIntensity(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingIntensity.CompareTo(right.MatchingIntensity);
        }
        public static int CompareMaxQuantScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MaxQuantScore().CompareTo(right.MaxQuantScore());
        }
        public static int CompareProductScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.ProductScore.CompareTo(right.ProductScore);
        }
        public static int ComparePrecursorScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.PrecursorScore.CompareTo(right.PrecursorScore);
        }
        public static int CompareMatchingProductsFraction(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingProductsFraction.CompareTo(right.MatchingProductsFraction);
        }
        public static int CompareMatchingProducts(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.MatchingWeightedProducts.CompareTo(right.MatchingWeightedProducts);
        }
        public static int CompareProteinScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.ProteinScore.CompareTo(right.ProteinScore);
        }
        public static int ComparePeptideScore(PeptideSpectrumMatch left, PeptideSpectrumMatch right)
        {
            return -left.PeptideScore.CompareTo(right.PeptideScore);
        }
    }
}
