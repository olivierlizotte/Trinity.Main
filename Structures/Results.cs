/*
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

        public void WriteInfoToCsv(bool light = false)
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
            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + "FragmentStats_" + (target ? "Targets" : "Decoy") + ".csv");
            writer.AddLine("  === Fragmentation of " + (target ? "Targets" : "Decoys") + " ===");
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
                writer.AddLine("    " + fragment.Name + ", Number of fragments = , " + nbFrag + ",   Intensity = ," + cumulIntensity + ", fragment matched [" + strPos.Substring(1) + "]");
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
                writer.AddLine("    " + fragment + ", Number of fragments = ," + nbFrag + ",   Intensity = ," + cumulIntensity);
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
                writer.AddLine("    " + fragment + ", Number of fragments = ," + nbFrag + ",   Intensity = ," + cumulIntensity);
            }
            writer.writeToFile();
        }

        public void Save()
        {
            GraphML graph = new GraphML();
            graph.Export(dbOptions.OutputFolder + "State.GraphML", this);
        }

        public void ExportFragments(PeptideSpectrumMatch psm)
        {
            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + psm.Peptide.Sequence + "_" + vsCSV.GetFileName_NoExtension(psm.Query.sample.sSDF) + "_" + psm.Query.precursor.Track.RT + ".csv");
            List<string> fragments = new List<string>();
            foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(psm.Peptide, psm.Query.precursor.Charge, dbOptions))
                {
                    if (fragment == match.fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }

            string title = "Theoretical Fragments";
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (string fragment in fragments)
                    title += "," + fragment + " ^" + charge;
            writer.AddLine(title);

            for (int i = 1; i <= psm.Peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(psm.Peptide, psm.Query.precursor.Charge, dbOptions))
                        {
                            if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.theoMz;
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                            line += ",";
                    }
                }
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (string fragment in fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(psm.Peptide, psm.Query.precursor.Charge, dbOptions))
                        {
                            if (fragment == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.theoMz;
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                            line += ",";
                    }
                }
                writer.AddLine(line);
            }

            title = "Observed Fragments Intensities";
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (string fragment in fragments)
                    title += "," + fragment + " ^" + charge; 
            writer.AddLine(title);

            for(int i = 1; i <= psm.Peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.obsIntensity;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (string fragment in fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.obsIntensity;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                writer.AddLine(line);
            }

            title = "Observed Fragments Mz";
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (string fragment in fragments)
                    title += "," + fragment + " ^" + charge;
            writer.AddLine(title);

            for (int i = 1; i <= psm.Peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.obsMz;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (string fragment in fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.obsMz;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                writer.AddLine(line);
            }

            title = "Error on Fragments";
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                foreach (string fragment in fragments)
                    title += "," + fragment + " ^" + charge;
            writer.AddLine(title);

            for (int i = 1; i <= psm.Peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.mass_diff;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                for (int charge = 1; charge <= psm.Query.precursor.Charge; charge++)
                {
                    foreach (string fragment in fragments)
                    {
                        bool found = false;
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                line += "," + match.mass_diff;
                                found = true;
                                break;
                            }
                        if (!found)
                            line += ",";
                    }
                }
                writer.AddLine(line);
            }
            writer.writeToFile();
        }

        public static double ComputePrecursorArea(List<PeptideSpectrumMatch> psms)
        {
            double fragSpectrumArea = 0;
            double lastTimeStamp = 0;
            foreach (PeptideSpectrumMatch psm in psms)
            {
                if (psm.Query.spectrum.PrecursorIntensity > 0 && lastTimeStamp > 0)
                    fragSpectrumArea += psm.Query.spectrum.PrecursorIntensity * (psm.Query.spectrum.RetentionTimeInMin - lastTimeStamp);
                lastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
            }
            return fragSpectrumArea;
        }

        public static Dictionary<PeptideSpectrumMatch, double> ComputeMsMsFactor(List<PeptideSpectrumMatch> psms)
        {
            Dictionary<PeptideSpectrumMatch, double> fragRatio = new Dictionary<PeptideSpectrumMatch, double>();
            double lastIntensity = 0;
            foreach (PeptideSpectrumMatch psm in psms)
            {
                if (psm.Query.spectrum.PrecursorIntensity > 0)
                {
                    double intensityFactor = 0;
                    if (psm.Query.spectrum.InjectionTime >= 119.999997317791)
                        lastIntensity = psm.Query.spectrum.PrecursorIntensity;
                    else
                        intensityFactor = (psm.Query.spectrum.PrecursorIntensity - lastIntensity) / lastIntensity;// psm.Query.spectrum.PrecursorIntensity;
                    fragRatio.Add(psm, intensityFactor);
                }
            }
            return fragRatio;
        }

        public List<ProductMatch> GetCommonSpectrum(List<PeptideSpectrumMatch> psms, Peptide peptide, int psmCharge, int nbProductsToKeep)
        {
            List<string> fragments = new List<string>();
            foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                int nbFrag = 0;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(peptide, psmCharge, dbOptions))
                {
                    nbFrag++;
                    if (fragment == match.fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }

            Dictionary<PeptideSpectrumMatch, double> MsMsFactor = ComputeMsMsFactor(psms);
//            double PrecursorArea = ComputePrecursorArea(psms);

            int nbExpectedPSM = 0;
            foreach (PeptideSpectrumMatch psm in psms)
                if (psm.MatchingProducts > 3)
                    nbExpectedPSM++;
            int nbCumuledFrag = 0;
            //Dictionary<ProductMatch, Dictionary<PeptideSpectrumMatch, int>> dicOfMatchToPSM = new Dictionary<ProductMatch, Dictionary<PeptideSpectrumMatch, int>>();
            List<ProductMatch> products = new List<ProductMatch>();
            for (int i = 1; i <= peptide.Length; i++)
            {
                for (int charge = 1; charge <= psmCharge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        int nbTimesSeen = 0;
                        ProductMatch pm = null;

                        double averageFragIntensity = 0;
                        //List<double> fragIntensities = new List<double>();
                        //double lastTimeStamp = 0;
                        foreach (PeptideSpectrumMatch psm in psms)
                        {
                            //if (psm.ProbabilityScore() > 0.02)
                            {
                                foreach (ProductMatch match in psm.AllProductMatches)
                                {
                                    if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                                    {
                                        pm = match;
                                        if (psm.Query.spectrum.PrecursorIntensity > 0)// && lastTimeStamp > 0)
                                        {
                                            averageFragIntensity += (match.obsIntensity + match.obsIntensity * MsMsFactor[psm]);
                                            //fragIntensities.Add(fragIntensity);
                                            //fragSpectrumArea.Add(psm.Query.spectrum.PrecursorIntensity * (psm.Query.spectrum.RetentionTimeInMin - lastTimeStamp));
                                            

                                            //cumul += fragIntensity;// / psm.Query.spectrum.PrecursorIntensity;//TODO Add precursor intensity                                            
                                            nbTimesSeen++;
               //                             if(!matchedPSMs.ContainsKey(psm))
                 //                               matchedPSMs.Add(psm, 0);
                                        }
                                        else
                                            Console.WriteLine("Null Intensity");
                                        //lastTimeStamp = psm.Query.spectrum.RetentionTimeInMin;
                                    }
                                }
                            }
                        }

                        if (pm != null)///* && cumul > 0*/ && fragIntensities.Count > 0)
                        {
                            ProductMatch savedPm = new ProductMatch(pm);
                            /*
                            double averageFragInt = 0;
                            double ratio = (fragSpectrumIntensities[k] * spectrumElapsedTime[k]) / precursorMaxIntensity;
                            double cumulRatio = 0;
                            for (int k = 0; k < fragIntensities.Count; k++)
                            {
                                if (ratio > 0.1)
                                {
                                    averageFragInt += ratio * fragIntensities[k];
                                    cumulRatio += ratio;
                                }
                            }
                            averageFragInt /= cumulRatio;//*/
                            savedPm.obsIntensity = averageFragIntensity;// cumul / (double)nbTimesSeen;
                            savedPm.weight = nbTimesSeen * savedPm.obsIntensity;// pm.obsIntensity;
                            products.Add(savedPm);
             //               dicOfMatchToPSM.Add(savedPm, matchedPSMs);
                            nbCumuledFrag++;
                        }/*
                        else
                        {
                            ProductMatch savedPm = new ProductMatch();
                            savedPm.charge = charge;
                            savedPm.fragment = fragment.Name;
                            savedPm.fragmentPos = i;
                            savedPm.mass_diff = 0;
                            savedPm.obsIntensity = 0;
                            foreach (ProductMatch matchTheo in dbOptions.fragments.ComputeFragments(peptide, 2, dbOptions))                                
                                    if (fragment.Name == matchTheo.fragment && matchTheo.fragmentPos == i && matchTheo.charge == charge)
                                    {
                                        savedPm.obsMz = matchTheo.theoMz;
                                        savedPm.theoMz = matchTheo.theoMz;
                                    }
                            savedPm.weight = 0;
                            products.Add(savedPm);
                            nbCumuledFrag++;
                        }//*/
                    }
                }
            }

            products.Sort(ProductMatch.AscendingWeightComparison);
            if (products.Count > nbProductsToKeep)
                products.RemoveRange(0, products.Count - nbProductsToKeep);//*/
            return products;
        }

        public List<double> GetCommonFragments(List<PeptideSpectrumMatch> psms, Peptide peptide, int psmCharge)
        {
            List<string> fragments = new List<string>();
            foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(peptide, psmCharge, dbOptions))
                {
                    if (fragment == match.fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }

            double cumulIntensity = 0;
            foreach (PeptideSpectrumMatch psm in psms)
                foreach (ProductMatch match in psm.AllProductMatches)
                    cumulIntensity += match.obsIntensity;
            
            int nbExpectedPSM = 0;
           foreach (PeptideSpectrumMatch psm in psms)
               if(psm.MatchingProducts > 3)
                   nbExpectedPSM ++;
           int nbCumuledFrag = 0;
            List<double> products = new List<double>();
            for (int i = 1; i <= peptide.Length; i++)
            {
                for (int charge = 1; charge <= psmCharge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        double cumul = 0.0;
                        int nbTimesSeen = 0;

                        foreach (PeptideSpectrumMatch psm in psms)
                        {
                            if(psm.MatchingProducts > 3)
                            {
                                foreach (ProductMatch match in psm.AllProductMatches)
                                {
                                    if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                                    {
                                        if (psm.Query.spectrum.PrecursorIntensity > 0)
                                        {
                                            cumul += match.obsIntensity / psm.Query.spectrum.PrecursorIntensity;//TODO Add precursor intensity
                                            nbTimesSeen++;
                                        }
                                        else
                                            Console.WriteLine("Null Intensity");
                                    }
                                }
                            }
                        }

                        if (cumul > 0 && nbTimesSeen >= nbExpectedPSM * 0.95)
                        {
                            products.Add(cumul / (double)nbTimesSeen);
                            nbCumuledFrag++;
                        }
                        else
                            products.Add(0);
                    }
                }
            }
            return products;
        }

        public List<double> GetFragments(PeptideSpectrumMatch psm, Peptide peptide, int psmCharge)
        {
            List<string> fragments = new List<string>();
            foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(peptide, psmCharge, dbOptions))
                {
                    if (fragment == match.fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }
            
            List<double> products = new List<double>();
            for (int i = 1; i <= peptide.Length; i++)
            {
                for (int charge = 1; charge <= psmCharge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        double cumul = 0.0;
                        int nbTimesSeen = 0;
                        
                        foreach (ProductMatch match in psm.AllProductMatches)
                            if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                            {
                                cumul += match.obsIntensity;
                                nbTimesSeen++;
                            }
                        products.Add(cumul);
                    }
                }
            }
            return products;
        }

        public List<double> GetFragments(List<PeptideSpectrumMatch> psms, Peptide peptide, int psmCharge, bool average, bool internalFragments, double minRatioIntensity)
        {
            List<string> fragments = new List<string>();
            foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(peptide, psmCharge, dbOptions))
                {
                    if (fragment == match.fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }

            double cumulIntensity = 0;
            foreach (PeptideSpectrumMatch psm in psms)
                foreach (ProductMatch match in psm.AllProductMatches)
                    cumulIntensity += match.obsIntensity;

            List<double> products = new List<double>();
            for (int i = 1; i <= peptide.Length; i++)
            {
                for (int charge = 1; charge <= psmCharge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        double cumul = 0.0;
                        int nbTimesSeen = 0;

                        foreach (PeptideSpectrumMatch psm in psms)
                            foreach (ProductMatch match in psm.AllProductMatches)
                                if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                                {
                                    cumul += match.obsIntensity;
                                    nbTimesSeen++;
                                }
                        if (cumul < cumulIntensity * minRatioIntensity)
                            cumul = 0;
                        if(average && nbTimesSeen > 0)
                            products.Add(cumul / (double)nbTimesSeen);
                        else
                            products.Add(cumul);
                    }
                    if (internalFragments)
                    {
                        foreach (string fragment in fragments)
                        {
                            double cumul = 0.0;
                            int nbTimesSeen = 0;
                            foreach (PeptideSpectrumMatch psm in psms)
                                foreach (ProductMatch match in psm.AllProductMatches)
                                    if (fragment == match.fragment && match.fragmentPos == i && match.charge == charge)
                                    {
                                        cumul += match.obsIntensity;
                                        nbTimesSeen++;
                                    }
                            if (average && nbTimesSeen > 0)
                                products.Add(cumul / (double)nbTimesSeen);
                            else
                                products.Add(cumul);
                        }
                    }
                }
            }
            return products;
        }

        public void ExportFragmentIntensities(List<PeptideSpectrumMatch> psms, Peptide peptide, int psmCharge, string fileName)
        {
            vsCSVWriter writer = new vsCSVWriter(fileName);
            List<string> fragments = new List<string>();
            foreach (string fragment in FragmentDictionary.Fragments.Keys)
            {
                bool found = false;
                foreach (ProductMatch match in dbOptions.fragments.ComputeFragments(peptide, psmCharge, dbOptions))
                {
                    if (fragment == match.fragment)
                    {
                        found = true;
                        break;
                    }
                }
                if (found)
                    fragments.Add(fragment);
            }

            string title = "Cumulated Product Intensities";
            for (int charge = 1; charge <= psmCharge; charge++)
                foreach (FragmentClass fragment in dbOptions.fragments)
                    title += "," + fragment.Name + " ^" + charge;
            for (int charge = 1; charge <= psmCharge; charge++)
                foreach (string fragment in fragments)
                    title += "," + fragment + " ^" + charge;
            writer.AddLine(title);

            for (int i = 1; i <= peptide.Length; i++)
            {
                string line = i.ToString();
                for (int charge = 1; charge <= psmCharge; charge++)
                {
                    foreach (FragmentClass fragment in dbOptions.fragments)
                    {
                        double cumul = 0.0;
                        foreach(PeptideSpectrumMatch psm in psms)
                            foreach (ProductMatch match in psm.AllProductMatches)
                                if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                                    cumul += match.obsIntensity;
                        line += "," + cumul;
                    }
                }
                for (int charge = 1; charge <= psmCharge; charge++)
                {
                    foreach (string fragment in fragments)
                    {
                        double cumul = 0.0;
                        foreach (PeptideSpectrumMatch psm in psms)
                            foreach (ProductMatch match in psm.AllProductMatches)
                                if (fragment == match.fragment && match.fragmentPos == i && match.charge == charge)
                                    cumul += match.obsIntensity;
                        line += "," + cumul;
                    }
                }
                writer.AddLine(line);
            }
            writer.writeToFile();
        }

        public void ExportFragmentIntensitiesForAllPSM(List<PeptideSpectrumMatch> psms, Peptide peptide, int psmCharge, string fileName)
        {
            vsCSVWriter writer = new vsCSVWriter(fileName);
            string title = "Retention Time";
            for (int i = 1; i <= peptide.Length; i++)
                for (int charge = 1; charge <= psmCharge; charge++)
                    foreach (FragmentClass fragment in dbOptions.fragments)
                        title += "," + i + fragment.Name + " ^" + charge;
            writer.AddLine(title);

            foreach (PeptideSpectrumMatch psm in psms)
            {
                string line = psm.Query.spectrum.RetentionTimeInMin.ToString();
                for (int i = 1; i <= peptide.Length; i++)
                {
                    for (int charge = 1; charge <= psmCharge; charge++)
                    {
                        foreach (FragmentClass fragment in dbOptions.fragments)
                        {
                            double cumul = 0.0;
                            foreach (ProductMatch match in psm.AllProductMatches)
                                if (fragment.Name == match.fragment && match.fragmentPos == i && match.charge == charge)
                                    cumul += match.obsIntensity;
                            line += "," + cumul;
                        }
                    }
                }
                writer.AddLine(line);
            }
            writer.writeToFile();
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
