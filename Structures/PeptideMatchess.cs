﻿/*
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
        public PeptideMatches() { }
        public PeptideMatches(PeptideMatch[] list) : base(list) { }

        public List<PeptideMatch> ComputeAtFDR(double desired_fdr)
        {
            if (uptimizer == null)
            {
                List<Comparison<PeptideMatch>> sorts = new List<Comparison<PeptideMatch>>();
                //sorts.Add(CompareMatchingProducts);
                //sorts.Add(CompareMatchingProductsFraction);
                //sorts.Add(CompareMatchingIntensityFraction);
                sorts.Add(ComparePeptideScore);
                sorts.Add(CompareProteinScore);
                //sorts.Add(CompareProductScore);
                sorts.Add(CompareCumulPrecursorScore);
                sorts.Add(CompareCumulPrecursorOptimizedScore);
                sorts.Add(CompareBestPrecursorScore);
                sorts.Add(CompareBestPrecursorOptimizedScore);
                sorts.Add(CompareScore);
                //sorts.Add(CompareNbCluster);
                //sorts.Add(ComparePrecursorMassError);
                uptimizer = new FDRizer<PeptideMatch>(this, sorts, null);
            }
            else
                uptimizer.ReStart();

            List<PeptideMatch> fdrList = uptimizer.Launch(desired_fdr);
            List<PeptideMatch> sortedProbability = new List<PeptideMatch>(this);
            sortedProbability.Sort(PeptideSearcher.DescendingOptimizedScoreComparison);
            sortedProbability = FDRizer<PeptideMatch>.ComputeAtFDR(sortedProbability, desired_fdr);
            if (sortedProbability.Count > fdrList.Count)
                fdrList = sortedProbability;

            return fdrList;
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
        public static int ComparePeptideScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).PeptideScore.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).PeptideScore);
        }
        public static int CompareProteinScore(PeptideMatch left, PeptideMatch right)
        {
            return -left.BestPrecursor().OptimizedBestPsm(left.peptide).ProteinScore.CompareTo(right.BestPrecursor().OptimizedBestPsm(right.peptide).ProteinScore);
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
            return -left.ScoreFct().CompareTo(right.ScoreFct());
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
