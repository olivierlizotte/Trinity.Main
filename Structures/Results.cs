﻿/*
 * Copyright 2013 Olivier Caron-Lizotte
 * olivierlizotte@gmail.com
 * Licensed under the MIT license: <http://www.opensource.org/licenses/mit-license.php>
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using Proteomics.Utilities;

namespace Trinity
{
    public class Result : GraphML_Node
    {
        public DBOptions dbOptions;
        public Samples samples;
        public Queries queries;
        public Precursors precursors;
        public Precursors matchedPrecursors { get; private set; }
        public GraphML_List<Cluster> clusters;
        public PeptideMatches peptides;
        public PeptideMatches peptideSequences;
        public ProteinGroupMatches proteins;

        public Result()
        {
            queries = new Queries();
            precursors = new Precursors();
            matchedPrecursors = new Precursors();
            clusters = new GraphML_List<Cluster>();
            peptides = new PeptideMatches(new PeptideMatch[0]);
            peptideSequences = new PeptideMatches(new PeptideMatch[0]);
            proteins = new ProteinGroupMatches();
            dbOptions = new DBOptions("");
            samples = new Samples();
        }

        public Result(Precursors precursors, GraphML_List<Cluster> clusters, PeptideMatches peptides, PeptideMatches peptideSequences, ProteinGroupMatches proteins, Queries queries, DBOptions dbOptions)
        {
            this.queries = queries;
            this.precursors = precursors;
            this.SetPrecursors(precursors);
            this.clusters = clusters;
            this.peptides = peptides;
            this.peptideSequences = peptideSequences;
            this.proteins = proteins;
            this.dbOptions = dbOptions;
        }

        public static Result Import(string filename)
        {
            GraphML g = new GraphML();
            return g.ImportGeneric(filename);
        }
        
        public long SetPrecursors(Precursors precursors)
        {
            this.precursors = precursors;
            this.matchedPrecursors = new Precursors();
            long nbTargets = 0;
            foreach (Precursor precursor in precursors)
                if (precursor.psms != null && precursor.psms.Count > 0)
                {
                    if (precursor.Target)
                        nbTargets++;
                    matchedPrecursors.Add(precursor);
                }
            return nbTargets;
        }
        /*
        public List<PeptideSpectrumMatch> GetPSMs(double fdr, double decoyOverTargetRatio = 1, int nbMaxPsm = 1)
        {
            return FDR.PSMs(clusters, fdr, decoyOverTargetRatio, nbMaxPsm);
        }//*/
        /*
        public List<PeptideMatch> GetPeptides(double fdr)
        {
            if (peptides != null)
                return FDR.Peptides(peptides, fdr);
            else
                return null;
        }

        public List<ProteinGroupMatch> GetProteins(double fdr)
        {
            if (proteins != null)
                return FDR.Proteins(proteins, fdr);
            else
                return null;
        }//*/
        
        public static int DescendingProteinScoreComparison(Precursor left, Precursor right)
        {
            PeptideSpectrumMatch psmLeft = left.OptimizedBestPsm();
            if (psmLeft != null)
            {
                PeptideSpectrumMatch psmRight = right.OptimizedBestPsm();
                if (psmRight != null)
                    return -(psmLeft.ProteinScore.CompareTo(psmRight.ProteinScore));
                else
                    return -1;
            }
            else return 1;
        }

        public static int DescendingIntensityFractionComparison(Precursor left, Precursor right)
        {
            return -(left.OptimizedBestPsm().MatchingIntensityFraction.CompareTo(right.OptimizedBestPsm().MatchingIntensityFraction));
        }

        public static int CountTargets(List<ITargetDecoy> elems)
        {
            int nbTarget = 0;
            foreach (ITargetDecoy t in elems)
                if (t.Target)
                    nbTarget++;
            return nbTarget;
        }

        public void WriteInfoToConsole(bool light = false)
        {
            int target = 0;
            foreach (Precursor precursor in matchedPrecursors)
                if (precursor.Target)
                    target++;
            Console.WriteLine("  ---  Number of precursors          : " + target + " targets [" + matchedPrecursors.Count + "]" + "  ---  ");

            target = 0;
            foreach (PeptideMatch peptide in peptides)
                if (peptide.Target)
                    target++;
            Console.WriteLine("  ---  Number of peptides            : " + target + " targets [" + peptides.Count + "]" + "  ---  ");

            target = 0;
            foreach (PeptideMatch peptide in peptideSequences)
                if (peptide.Target)
                    target++;
            Console.WriteLine("  ---  Number of peptide sequences   : " + target + " targets [" + peptideSequences.Count + "]" + "  ---  ");

            target = 0;
            foreach (ProteinGroupMatch protein in proteins)
                if (protein.Target)
                    target++;
            Console.WriteLine("  ---  Number of proteins            : " + target + " targets [" + proteins.Count + "]" + "  ---  ");
            if (!light)
            {
                WriteFragmentation(true);
                WriteFragmentation(false);
            }
        }

        public void WriteFragmentation(bool target)
        {
            Console.WriteLine("  === Fragmentation of " + (target ? "Targets" : "Decoys") + " ===");
            foreach (FragmentClass fragment in dbOptions.fragments)
            {
                double cumulIntensity = 0;
                int nbFrag = 0;
                Dictionary<int, int> positions = new Dictionary<int, int>();
                foreach (Precursor precursor in matchedPrecursors)
                {
                    PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                    if (psm.Target == target)
                    {
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment.Name == match.fragment)
                            {
                                nbFrag++;
                                if (!positions.ContainsKey(match.fragmentPos))
                                    positions.Add(match.fragmentPos, 1);
                                else
                                    positions[match.fragmentPos]++;
                                cumulIntensity += match.obsIntensity;
                            }
                    }
                }
                string strPos = "";
                if (positions.Count > 0)
                    foreach (int key in positions.Keys)
                        strPos += "|" + key + ":" + positions[key];
                else
                    strPos += ",";
                Console.WriteLine("    " + fragment.Name + ", Number of fragments = , " + nbFrag + ",   Intensity = ," + cumulIntensity);
                Console.WriteLine("     fragment matched [" + strPos.Substring(1) + "]");
            }
            foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                double cumulIntensity = 0;
                int nbFrag = 0;
                foreach (Precursor precursor in matchedPrecursors)
                {
                    PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                    if (psm.Target == target)
                    {
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.fragment)
                            {
                                nbFrag++;
                                cumulIntensity += match.obsIntensity;
                            }
                    }
                }
                Console.WriteLine("    " + fragment + ", Number of fragments = ," + nbFrag + ",   Intensity = ," + cumulIntensity);
            }
            foreach (string fragment in FragmentDictionary.AAFragments.Keys)
            {
                double cumulIntensity = 0;
                int nbFrag = 0;
                foreach (Precursor precursor in matchedPrecursors)
                {
                    PeptideSpectrumMatch psm = precursor.OptimizedBestPsm();
                    if (psm.Target == target)
                    {
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.fragment)
                            {
                                nbFrag++;
                                cumulIntensity += match.obsIntensity;
                            }
                    }
                }
                Console.WriteLine("    " + fragment + ", Number of fragments = ," + nbFrag + ",   Intensity = ," + cumulIntensity);
            }
        }

        public void Save()
        {
            GraphML graph = new GraphML();
            graph.Export(dbOptions.OutputFolder + "State.GraphML", this);
        }

        public void Export(double fdr, string keyword = "", bool onlyPrecursors = false)
        {
            Console.WriteLine("Exporting at " + (fdr * 100) + "% FDR (Decoy/Target)...");

            List<Precursor> prec = null;
            if (precursors != null)
            {
                if (matchedPrecursors == null)
                {
                    this.matchedPrecursors = new Precursors();
                    foreach (Precursor precursor in precursors)
                        if (precursor.psms.Count > 0)
                            matchedPrecursors.Add(precursor);
                }/*
                List<Precursor> prec = FDR.PrecursorsV2(precursors, fdr, 1);
                Sol.CONSOLE.OutputLine(">   " + prec.Count + " Precursors");
                MSSearcher.Export(dbOptions.outputFolder + keyword + "precursors.csv", prec);//*/
                /*
                prec = Optimizer.PrecursorOptimizer(matchedPrecursors, fdr);
                Sol.CONSOLE.OutputLine(">   " + prec.Count + " Optimized Precursors");
                MSSearcher.Export(dbOptions.outputFolder + keyword + "Optimized_precursors.csv", prec);//*/

                prec = matchedPrecursors.ComputeAtFDR(fdr);
                Console.WriteLine(">   " + prec.Count + " Uptimized V5 Precursors");
                MSSearcher.Export(dbOptions.OutputFolder + keyword + "UptimizedV5_precursors.csv", prec);
            }

            if (!onlyPrecursors)
            {
                if (queries != null)
                {
                    List<Query> qs = queries.ComputeAtFDR(fdr);
                    Console.WriteLine(">   " + qs.Count + " PSMs (Top 10)");
                    MSSearcher.Export(dbOptions.OutputFolder + keyword + "queries.csv", qs);
                }
                if (clusters != null)
                {
                    //List<PeptideSpectrumMatch> psms = FDR.PSMs(clusters, fdr, 1, 10);
                    //Console.WriteLine(">   " + psms.Count + " PSMs (Top 10)");
                    //MSSearcher.Export(dbOptions.outputFolder + keyword + "psms_Top10.csv", psms);

                    //psms = FDR.PSMs(clusters, fdr, 1, 1);
                    //Console.WriteLine(">   " + psms.Count + " PSMs");
                    //MSSearcher.Export(dbOptions.outputFolder + keyword + "psms_Best.csv", psms);
                }

                if (peptides != null)
                {
                    List<PeptideMatch> pep = peptideSequences.ComputeAtFDR(fdr);
                    Console.WriteLine(">   " + pep.Count + " Peptides Sequences (Version 5)");
                    PeptideSearcher.Export(dbOptions.OutputFolder + keyword + "peptideSequencesV5_.csv", pep);

                    PeptideSearcher sr = new PeptideSearcher(dbOptions);
                    PeptideMatches seqs = sr.Search(clusters, prec, false);
                    pep = seqs.ComputeAtFDR(fdr);
                    Console.WriteLine(">   " + pep.Count + " Peptides Sequences (Version 5b)");
                    PeptideSearcher.Export(dbOptions.OutputFolder + keyword + "peptideSequencesV5b_PrecursorFDRed.csv", pep);

                    pep = peptides.ComputeAtFDR(fdr);
                    Console.WriteLine(">   " + pep.Count + " Peptides (Version 5)");
                    PeptideSearcher.Export(dbOptions.OutputFolder + keyword + "peptidesV5_.csv", pep);

                }
                if (proteins != null)
                {
                    List<ProteinGroupMatch> prots = proteins.ComputeAtFDR(fdr);
                    Console.WriteLine(">   " + prots.Count + " Proteins");
                    ProteinSearcher.Export(dbOptions.OutputFolder + keyword + "proteins_.csv", prots);
                }
            }
        }
    }    
}