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
    public class PeptideMatches : GraphML_List<PeptideMatch>
    {
        private FDRizer<PeptideMatch> uptimizer;
        private FDRizer2<PeptideMatch> uptimizer2;

        public PeptideMatches() { }
        public PeptideMatches(IEnumerable<PeptideMatch> list) : base(list) { }

        public List<PeptideMatch> ComputeAtFDR(double desired_fdr)
        {
            if (uptimizer == null)
            {
                List<Comparison<PeptideMatch>> sorts = new List<Comparison<PeptideMatch>>();
                //sorts.Add(CompareMatchingProducts);
                //sorts.Add(CompareMatchingProductsFraction);
                //sorts.Add(CompareMatchingIntensityFraction);
         
                //sorts.Add(CompareProductScore);
                //                sorts.Add(CompareCumulPrecursorScore);
                sorts.Add(CompareCumulPrecursorOptimizedScore);
                //                sorts.Add(CompareBestPrecursorScore);
                sorts.Add(CompareBestPrecursorOptimizedScore);
                sorts.Add(CompareScore);
                //sorts.Add(CompareNbCluster);
                //sorts.Add(ComparePrecursorMassError);
                uptimizer = new FDRizer<PeptideMatch>(this, sorts, null);
            }
            else
                uptimizer.ReStart();

            if (uptimizer2 == null)
            {
                List<Comparison<PeptideMatch>> sorts = new List<Comparison<PeptideMatch>>();
                //sorts.Add(CompareMatchingProducts);
                //sorts.Add(CompareMatchingProductsFraction);
                //sorts.Add(CompareMatchingIntensityFraction);
         
                //sorts.Add(CompareProductScore);
//                sorts.Add(CompareCumulPrecursorScore);
                sorts.Add(CompareCumulPrecursorOptimizedScore);
//                sorts.Add(CompareBestPrecursorScore);
                sorts.Add(CompareBestPrecursorOptimizedScore);
                sorts.Add(CompareScore);
                //sorts.Add(CompareNbCluster);
                //sorts.Add(ComparePrecursorMassError);
                uptimizer2 = new FDRizer2<PeptideMatch>(this, sorts, null);
            }
            else
                uptimizer2.ReStart();

            List<PeptideMatch> fdrList = uptimizer.Launch(desired_fdr);
            List<PeptideMatch> fdrList2 = uptimizer2.Launch(desired_fdr);

            if (fdrList.Count < fdrList2.Count)
                return fdrList2;
            else
                return fdrList;
            //List<PeptideMatch> sortedProbability = new List<PeptideMatch>(this);
            //sortedProbability.Sort(PeptideSearcher.DescendingOptimizedScoreComparison);
            //sortedProbability = FDRizer<PeptideMatch>.ComputeAtFDR(sortedProbability, desired_fdr);
            //if (sortedProbability.Count > fdrList.Count)
            //    fdrList = sortedProbability;
            //return fdrList;
        }

        public static int CompareMatchingProducts(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).MatchingProducts.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).MatchingProducts);
        }
        public static int CompareMatchingProductsFraction(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).MatchingProductsFraction.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).MatchingProductsFraction);
        }
        public static int CompareMatchingIntensityFraction(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).MatchingIntensityFraction.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).MatchingIntensityFraction);
        }
        public static int CompareProductScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).ProductScore.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).ProductScore);
        }
        public static int CompareCumulPrecursorOptimizedScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.CumulPrecursorOptimizedScore().CompareTo(right.CumulPrecursorOptimizedScore());
        }
        public static int CompareCumulPrecursorScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.CumulPrecursorScore().CompareTo(right.CumulPrecursorScore());
        }
        public static int CompareBestPrecursorOptimizedScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursorOptimizedScore().CompareTo(right.BestPrecursorOptimizedScore());
        }
        public static int CompareBestPrecursorScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursorScore().CompareTo(right.BestPrecursorScore());
        }
        public static int CompareScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }
        public static int CompareNbCluster(PeptideMatch left, PeptideMatch right)
        {
            return -left.clusters.Count.CompareTo(right.clusters.Count);
        }
        public static int ComparePrecursorMassError(PeptideMatch left, PeptideMatch right)
        {
            return left.GetPrecursorMassError().CompareTo(right.GetPrecursorMassError());
        }
    }
}
