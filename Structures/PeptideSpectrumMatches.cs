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

        public void OptimizePSMScoreRatios(DBOptions options, double desired_fdr, Result results)
        {
            long bestNbTargets = 0;
            double bestProtein = 0.1;
            double bestPeptide = 0.2;
            double bestFragments = 0.4;
            double bestIntensities = 0.1;
            double bestPrecursor = 0.2;

            double incr = 0.05;
            for (options.dProtein = 0.1; options.dProtein < 0.3; options.dProtein += incr)
                for (options.dPeptideScore = 0.1; options.dPeptideScore < 0.4; options.dPeptideScore += incr)
                    for (options.dPrecursor = 0.1; options.dPrecursor < 0.5; options.dPrecursor += incr)
                        for (options.dMatchingProductFraction = 0.1; options.dMatchingProductFraction < 1 - options.dPrecursor - options.dPeptideScore - options.dProtein; options.dMatchingProductFraction += incr)
                        {
                            double cumul = options.dPrecursor + options.dPeptideScore + options.dProtein + options.dMatchingProductFraction;
                            for (options.dIntensityFraction = 0.1; options.dIntensityFraction < 1 - cumul; options.dIntensityFraction += incr)
                            {
                                if (cumul + options.dIntensityFraction == 1)
                                {
                                    long nbTargets = 0;
                                    foreach (Precursor precursor in results.matchedPrecursors)
                                        if (precursor.Target)
                                            nbTargets++;
                                    if (nbTargets > bestNbTargets)
                                    {
                                        bestNbTargets = nbTargets;
                                        bestProtein = options.dProtein;
                                        bestPeptide = options.dPeptideScore;
                                        bestPrecursor = options.dPrecursor;
                                        bestIntensities = options.dIntensityFraction;
                                        bestFragments = options.dMatchingProductFraction;
                                    }
                                }
                            }

                        }

            options.dProtein = bestProtein;
            options.dPeptideScore = bestPeptide;
            options.dPrecursor = bestPrecursor;
            options.dIntensityFraction = bestIntensities;
            options.dMatchingProductFraction = bestFragments;
            
            Console.WriteLine("New score ratios   ----------------------------------------------------------------------- ");
            Console.WriteLine("    PeptideSpectrumMatch.dPrecursor:                     " + options.dPrecursor);
            Console.WriteLine("    PeptideSpectrumMatch.dMatchingProductFraction:       " + options.dMatchingProductFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dIntensityFraction:             " + options.dIntensityFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dPeptideScore:                  " + options.dPeptideScore);
            Console.WriteLine("    PeptideSpectrumMatch.dProtein:                       " + options.dProtein);
            Console.WriteLine("------------------------------------------------------------------------------------------ ");
        }

        public void OptimizePSMScoreRatios_PREVIOUSVersion(DBOptions options, double desired_fdr)
        {
            //TODO find a Max flow approach to modelize the score optimization routine to maximize Targets versus Decoys
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
            options.dPrecursor += ratioPrecursorPrecision / totalRatios;
            options.dMatchingProductFraction += ratioFragments / totalRatios;
            options.dIntensityFraction += ratioIntensities / totalRatios;
            options.dPeptideScore += ratioPeptideScore / totalRatios;

            double somme = options.dPrecursor + options.dMatchingProductFraction + options.dIntensityFraction + options.dPeptideScore;
            options.dPrecursor /= somme;
            options.dMatchingProductFraction /= somme;
            options.dIntensityFraction /= somme;
            options.dPeptideScore /= somme;

            //options.dPrecursor                 = ratioPrecursorPrecision   / totalRatios;
            //options.dMatchingProductFraction   = ratioFragments            / totalRatios;
            //options.dIntensityFraction         = ratioIntensities          / totalRatios;
            //options.dPeptideScore              = ratioPeptideScore         / totalRatios;
            Console.WriteLine("New score ratios  [" + totalRatios + " total ratios] ------------------------------------- ");
            Console.WriteLine("    PeptideSpectrumMatch.dPrecursor:                     " + options.dPrecursor);
            Console.WriteLine("    PeptideSpectrumMatch.dMatchingProductFraction:       " + options.dMatchingProductFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dIntensityFraction:             " + options.dIntensityFraction);
            Console.WriteLine("    PeptideSpectrumMatch.dPeptideScore:                  " + options.dPeptideScore);
            Console.WriteLine("------------------------------------------------------------------------------------------ ");
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
