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
    /// <summary>
    /// Peptide sequence match to at least one spectrum
    /// Holds the list of clusters associated to this sequence
    /// </summary>
    public class PeptideMatch : GraphML_Node, ITargetDecoy
    {
        public bool Decoy
        { get { return peptide.Decoy; } }

        public bool Target
        { get { return !peptide.Decoy; } }


        public GraphML_List<Cluster> clusters;
        public Peptide peptide;
        public PeptideMatch()
        {
            clusters = new GraphML_List<Cluster>();
        }
        public PeptideMatch(Peptide peptide)
        {
            this.peptide = peptide;
            clusters = new GraphML_List<Cluster>();
        }
        public void AddOnlyOnce(Cluster cluster)
        {
            if (!clusters.Contains(cluster))
                clusters.Add(cluster);
        }

        public static int DescendingOptimizedScoreComparison(PeptideMatch left, PeptideMatch right)
        {
            return -(left.ProbabilityScore().CompareTo(right.ProbabilityScore()));
        }

        public static int DescendingPrecursorScoreComparison(PeptideMatch left, PeptideMatch right)
        {
            return -(left.CumulPrecursorScore().CompareTo(right.CumulPrecursorScore()));
        }

        public double ProbabilityScore()
        {
            double score = 0;
            foreach (Cluster cluster in clusters)
                score += (1 - score) * cluster.ProbabilityScore(peptide);
            return score;
        }


        public double BestPrecursorScore()
        {
            double score = 0;
            foreach(Cluster cluster in clusters)
                foreach(clCondition condition in cluster.conditions)
                    foreach(clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            if (precursor.ProbabilityScore(peptide) > score)
                                score = precursor.ProbabilityScore(peptide);
            return score;
        }
        
        public double BestPrecursorOptimizedScore()
        {
            double score = 0;
            foreach (Cluster cluster in clusters)
                foreach (clCondition condition in cluster.conditions)
                    foreach (clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            if (precursor.ProbabilityScore() > score)
                                score = precursor.ProbabilityScore(peptide);
            return score;
        }

        public Precursor BestPrecursor(bool checkMods = false)
        {
            double score = 0;
            Precursor best = null;
            foreach (Cluster cluster in clusters)
            {
                Precursor tmp = cluster.OptimizedBestPrecursor(peptide, checkMods);
                if (tmp != null)
                {
                    double tmpScore = tmp.ProbabilityScore(peptide, checkMods);
                    if (tmpScore > score)
                    {
                        score = tmpScore;
                        best = tmp;
                    }
                }
            }
            return best;
        }
        
        public double CumulPrecursorScore()
        {
            double score = 0;
            foreach (Cluster cluster in clusters)
                foreach (clCondition condition in cluster.conditions)
                    foreach (clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            score += precursor.ProbabilityScore(peptide);
            return score;
        }

        public double CumulPrecursorOptimizedScore()
        {
            double score = 0;
            foreach (Cluster cluster in clusters)
                foreach (clCondition condition in cluster.conditions)
                    foreach (clReplicate replicate in condition.replicates)
                        foreach (Precursor precursor in replicate.precursors)
                            score += precursor.ProbabilityScore(peptide);
            return score;
        }//*/

        public double GetPrecursorMassError()
        {
            double cumul = 0;
            int nbSeen = 0;
            foreach(Cluster cluster in clusters)
            {
                cumul += cluster.GetPrecursorMassError(peptide);
                nbSeen++;
            }
            return cumul / (double)nbSeen;
        }
    }

    /// <summary>
    /// Methods to parse the set of clusters and precursors to get a list of expressed peptide sequences
    /// </summary>
    public class PeptideSearcher
    {
        public DBOptions options;
        public PeptideSearcher(DBOptions options)
        {
            this.options = options;
        }
        
        public static int AscendingSpectrumNumberComparison(Cluster left, Cluster right)
        {
            return left.Rt().CompareTo(right.Rt());
        }
        
        public static int DescendingOptimizedScoreComparison(PeptideMatch left, PeptideMatch right)
        {
            return -left.ProbabilityScore().CompareTo(right.ProbabilityScore());
        }

        
        public PeptideMatches Search(List<Cluster> clusters, List<Precursor> precursors, bool DiffByMod = true)
        {
            Console.WriteLine("Creating the list of peptide found...");
            Dictionary<string, PeptideMatch> peptideMatches = new Dictionary<string, PeptideMatch>();            

            precursors.Sort(Precursor.CompareProbabilityScore);//.DescendingScoreComparison);
            foreach (Precursor precursor in precursors)
                foreach(PeptideSpectrumMatch psm in precursor.OptimizedBestPsms())
                {
                    string seq = (DiffByMod ? psm.Peptide.Sequence : psm.Peptide.BaseSequence);
                    if (!peptideMatches.ContainsKey(seq))
                        peptideMatches.Add(seq, new PeptideMatch(psm.Peptide));
                    else if (psm.Peptide.Target)
                        peptideMatches[seq].peptide = psm.Peptide;
                }

            foreach (Cluster cluster in clusters)
                foreach(PeptideSpectrumMatch psm in cluster.OptimizedBestPsms())
                {
                    string seq = (DiffByMod ? psm.Peptide.Sequence : psm.Peptide.BaseSequence);
                    if(peptideMatches.ContainsKey(seq))
                        peptideMatches[seq].AddOnlyOnce(cluster);
                }

            PeptideMatches TotalList = new PeptideMatches();
            foreach (PeptideMatch match in peptideMatches.Values)
                if (match.clusters.Count > 0)
                    TotalList.Add(match);
            //PeptideMatchess TotalList = new PeptideMatchess(peptideMatches.Values.ToArray<PeptideMatch>());//<PeptideMatch>());
            Console.WriteLine(TotalList.Count + " distinct peptides (based on sequence" + (DiffByMod ? " and modifications)" : ")"));
            return TotalList;
        }

        public static void Export(string filename, List<PeptideMatch> peptides)
        {
            vsCSVWriter writer = new vsCSVWriter(filename);
            writer.AddLine("Sequence,Variable Modification,Score,Decoy,Precursor Mass Error");
            foreach (PeptideMatch pm in peptides)
                writer.AddLine(pm.peptide.BaseSequence + "," + pm.peptide.Sequence + "," + pm.ProbabilityScore() + "," + pm.peptide.Decoy + "," + pm.GetPrecursorMassError());
            writer.writeToFile();
        }
    }
}
