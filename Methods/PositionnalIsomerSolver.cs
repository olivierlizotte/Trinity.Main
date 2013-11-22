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

namespace Trinity.UnitTest
{
    public class CharacterizedPrecursor : PrecursorIon
    {
        public Peptide Peptide;
        public PeptideSpectrumMatches Psms;
        public List<ProductMatch> Fragments;
        public CharacterizedPrecursor(Sample sample, DBOptions dbOptions, Peptide peptide, IEnumerable<Query> queries, double mz, int charge = -1):base(sample, queries, mz, charge)
        {
            this.Peptide = peptide;
            Psms = new PeptideSpectrumMatches();
            foreach(Query query in queries)
                if(query.sample == sample)
                    Psms.Add(new PeptideSpectrumMatch(query, peptide, dbOptions));
            Psms.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
            Fragments = Psms.GetCombinedSpectrum(dbOptions);
            Fragments.Sort(ProductMatch.DescendingWeightComparison);
        }
    }

    public class MixedPrecursor : PrecursorIon
    {
        public Dictionary<Sample, ElutionCurve> PeptideRatios;
        public MixedPrecursor(Sample sample, IEnumerable<Query> queries, double mz):base(sample, queries, mz, -1)
        {
        }
    }

    public class PrecursorIon
    {
        public ElutionCurve eCurve;
        public double MZ;
        public int Charge;
        public List<Query> Queries;
        public Sample Sample;
        public PrecursorIon(Sample sample, IEnumerable<Query> queries, double mz, int charge)
        {
            this.MZ = mz;
            this.Charge = charge;
            this.Queries = queries.ToList<Query>();
            this.Queries.Sort(Query.AscendingRetentionTimeComparison);            
            this.eCurve = ElutionCurve.Create(this.Queries);
            this.Sample = sample;
        }
    }

    public class PositionnalIsomerSolver
    {
        private static DBOptions CreateOptions(string fastaFile, string outputFolder, IConSol consol)
        {
            DBOptions dbOptions = new DBOptions(fastaFile, consol);
            dbOptions.precursorMassTolerance = new MassTolerance(8, MassToleranceUnits.ppm);
            dbOptions.productMassTolerance = new MassTolerance(20, MassToleranceUnits.ppm);
            dbOptions.MaximumPeptideMass = 200000;
            dbOptions.OutputFolder = outputFolder;

            ProteaseDictionary proteases = ProteaseDictionary.Instance;
            dbOptions.DigestionEnzyme = proteases["no enzyme"];
            dbOptions.NoEnzymeSearch = false;
            dbOptions.DecoyFusion = false;
            dbOptions.MaximumNumberOfFragmentsPerSpectrum = 400;
            dbOptions.ToleratedMissedCleavages = 200;
            dbOptions.MinimumPeptideLength = 5;
            dbOptions.MaximumPeptideLength = 300;

            GraphML_List<Modification> fixMods = new GraphML_List<Modification>();
            dbOptions.fixedModifications = fixMods;

            GraphML_List<Modification> varMods = new GraphML_List<Modification>();
            foreach (string strMod in ModificationDictionary.Instance.Keys)
                varMods.Add(ModificationDictionary.Instance[strMod]);

            dbOptions.maximumVariableModificationIsoforms = 1024;
            dbOptions.variableModifications = varMods;

            dbOptions.addFragmentLoss = false;
            dbOptions.addFragmentMods = false;
            dbOptions.fragments = new Fragments();

            dbOptions.fragments.Add(new FragmentA());
            dbOptions.fragments.Add(new FragmentB());
            dbOptions.fragments.Add(new FragmentC());
            dbOptions.fragments.Add(new FragmentX());
            dbOptions.fragments.Add(new FragmentY());
            dbOptions.fragments.Add(new FragmentZ());

            dbOptions.SaveMS1Peaks = true;
            dbOptions.SaveMSMSPeaks = true;
            dbOptions.LoadSpectraIfFound = true;

            dbOptions.NbPSMToKeep = 100;
            return dbOptions;
        }

        public static void Solve(string[] spikedRaws, string[] mixedRaws, string fastaFile, string outputFolder, IConSol conSol)
        {
            DBOptions dbOptions = CreateOptions(fastaFile, outputFolder, conSol);
            Samples SpikedSamples = new Samples(dbOptions);
            for (int i = 0; i < spikedRaws.Length; i++)
                SpikedSamples.Add(new Sample(i + 1, 1, 1, spikedRaws[i], spikedRaws[i], 0, ""));

            //Precompute Spiked peptide identifications
            Result spikedResult = Propheus.Start(dbOptions, SpikedSamples, false, false, true);
            
            Samples MixedSamples = new Samples(dbOptions);
            for (int i = 0; i < mixedRaws.Length; i++)
                MixedSamples.Add(new Sample(i + 1, 1, 1, mixedRaws[i], mixedRaws[i], 0, ""));

            //Precompute Mixed peptide identifications
            Result mixedResult = Propheus.Start(dbOptions, MixedSamples, false, false, true);

            //!!!! Compute spiked peptide profiles and list of fragments
            //!!!! Compute Normalization vectors
            Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> characterizedPeptides = GetSpikedPrecursors(SpikedSamples, spikedResult, dbOptions);
                        
            //Get the list of precursors to characterize
            foreach (Sample mixedSample in MixedSamples)
            {
                Dictionary<double, MixedPrecursor> mixedPrecursors = GetMixedPrecursors(mixedSample, mixedResult, dbOptions, characterizedPeptides);

                foreach (double keyMz in mixedPrecursors.Keys)
                {
                    // Compute Max Flow for this precursor
                    Dictionary<Sample, double> ratios = GetRatiosFromMaxFlow(dbOptions, characterizedPeptides, mixedPrecursors[keyMz]);
                    
                    // Export results in a file
                    vsCSVWriter writerCumul = new vsCSVWriter(dbOptions.OutputFolder + vsCSV.GetFileName_NoExtension(mixedSample.sSDF) + "_" + keyMz + "MZ.csv");
                    writerCumul.AddLine("Sequence,Area");
                    writerCumul.AddLine("Total," + mixedPrecursors[keyMz].eCurve.Area);
                    foreach (Sample spikedSample in ratios.Keys)
                        writerCumul.AddLine(characterizedPeptides[keyMz][spikedSample].Peptide.Sequence + "," + ratios[spikedSample]);
                    writerCumul.WriteToFile();
                }
            }
        }

        public static Dictionary<double, List<Query>> GetPrecursors(Result result, Sample sample, DBOptions dbOptions, IEnumerable<double> keys = null)
        {
            Dictionary<double, List<Query>> DicOfSpectrumMasses = new Dictionary<double, List<Query>>();
            foreach (Query query in result.queries)
            {
                if (query.sample == sample)
                {
                    double foundKey = query.spectrum.PrecursorMZ;
                    if(keys != null)
                    {
                        foreach (double key in keys)
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(query.spectrum.PrecursorMZ, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                foundKey = key;                            
                    }else if(query.psms.Count > 0)
                        foundKey = Numerics.MZFromMass(query.psms[0].Peptide.MonoisotopicMass, query.spectrum.PrecursorCharge);

                    if (!DicOfSpectrumMasses.ContainsKey(foundKey))
                    {
                        List<Query> listOfSpectrum = new List<Query>();
                        listOfSpectrum.Add(query);
                        DicOfSpectrumMasses.Add(foundKey, listOfSpectrum);
                    }
                    else
                        DicOfSpectrumMasses[foundKey].Add(query);
                }
            }
            return DicOfSpectrumMasses;
        }

        public static Dictionary<double, MixedPrecursor> GetMixedPrecursors(Sample mixedSample, Result mixedResult, DBOptions dbOptions, Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> charPeptides)
        {
            Dictionary<double, List<Query>> DicOfSpectrumMasses = GetPrecursors(mixedResult, mixedSample, dbOptions, charPeptides.Keys);
            Dictionary<double, MixedPrecursor> DicOfMixedPrecursor = new Dictionary<double,MixedPrecursor>();
            foreach (double key in DicOfSpectrumMasses.Keys)
            {
                if (charPeptides.ContainsKey(key))
                {
                    MixedPrecursor mixedPrecursor = new MixedPrecursor(mixedSample, DicOfSpectrumMasses[key], key);
                    DicOfMixedPrecursor.Add(key, mixedPrecursor);
                }
            }
            return DicOfMixedPrecursor;            
        }

        public static Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> GetSpikedPrecursors(Samples spikedSamples, Result spikedResult, DBOptions dbOptions)
        {
            Dictionary<double, Dictionary<Query, int>> mzKeys = new Dictionary<double, Dictionary<Query, int>>();
            foreach(Query query in spikedResult.queries)
            {
                foreach(PeptideSpectrumMatch psm in query.psms)
                {
                    double mz = Numerics.MZFromMass(psm.Peptide.MonoisotopicMass, query.spectrum.PrecursorCharge);
                    if(!mzKeys.ContainsKey(mz))
                    {
                        bool found = false;
                        foreach(double key in mzKeys.Keys)
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(mz, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                            {
                                mz = key;
                                found = true;
                            }
                        if(!found)
                            mzKeys.Add(mz, new Dictionary<Query, int>());
                    }
                    if(!mzKeys[mz].ContainsKey(query))
                        mzKeys[mz].Add(query, 1);
                }
            }
                
            Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> spikes = new Dictionary<double,Dictionary<Sample,CharacterizedPrecursor>>();
            foreach(Sample spikedSample in spikedSamples)
            {
                Dictionary<double, List<Query>> DicOfSpectrumMasses = GetPrecursors(spikedResult, spikedSample, dbOptions, mzKeys.Keys);                
                foreach(double mzKey in DicOfSpectrumMasses.Keys)
                {
                    if (mzKeys.ContainsKey(mzKey))
                    {
                        //Pick the best PSM for each sample/precursor pair
                        Dictionary<Peptide, double> DicOfProbabilityScores = new Dictionary<Peptide, double>();

                        foreach (Query query in mzKeys[mzKey].Keys)
                            if (query.sample == spikedSample)
                            {
                                foreach (PeptideSpectrumMatch psm in query.psms)
                                    if (!DicOfProbabilityScores.ContainsKey(psm.Peptide))
                                        DicOfProbabilityScores.Add(psm.Peptide, psm.ProbabilityScore());
                                    else
                                        DicOfProbabilityScores[psm.Peptide] += psm.ProbabilityScore();
                            }

                        Peptide bestPeptide = null;
                        double bestScore = double.MinValue;
                        foreach (Peptide keyPep in DicOfProbabilityScores.Keys)
                            if (DicOfProbabilityScores[keyPep] > bestScore)
                            {
                                bestScore = DicOfProbabilityScores[keyPep];
                                bestPeptide = keyPep;
                            }
                        if (bestPeptide != null)
                        {
                            CharacterizedPrecursor cPrec = new CharacterizedPrecursor(spikedSample, dbOptions, bestPeptide, mzKeys[mzKey].Keys, mzKey);
                            if (!spikes.ContainsKey(mzKey))
                                spikes.Add(mzKey, new Dictionary<Sample, CharacterizedPrecursor>());
                            if (!spikes[mzKey].ContainsKey(spikedSample))
                                spikes[mzKey].Add(spikedSample, cPrec);
                            else
                                Console.WriteLine("Twice??");
                        }
                    }
                }//End of foreach mzKey
            }//End of foreach spiked sample
            return spikes;
        }

        public static Dictionary<double, int> GetCommonFragmentMz(IEnumerable<CharacterizedPrecursor> precursors, int nbFragmentPerPrec)
        {
            Dictionary<double, int> DicOfFragments = new Dictionary<double,int>();
            foreach(CharacterizedPrecursor cPrec in precursors)
            {
                for(int i = 0; i < nbFragmentPerPrec && i < cPrec.Fragments.Count; i++)
                    if(!DicOfFragments.ContainsKey(cPrec.Fragments[i].theoMz))
                        DicOfFragments.Add(cPrec.Fragments[i].theoMz, 1);
                    else
                        DicOfFragments[cPrec.Fragments[i].theoMz] ++;
            }
            return DicOfFragments;
        }
        
        public static Dictionary<Sample, double> GetRatiosFromMaxFlow(DBOptions dbOptions, Dictionary<double, Dictionary<Sample, CharacterizedPrecursor>> spikes, MixedPrecursor mixedPrecursor)        
        {
            int minProducts = 5;
            int maxProducts = 14;
            Dictionary<Sample, double> values = new Dictionary<Sample,double>();
            
            for (int nbProductsToKeep = minProducts; nbProductsToKeep <= maxProducts; nbProductsToKeep++)
            {
                Dictionary<Sample, CharacterizedPrecursor> Isomers = new Dictionary<Sample, CharacterizedPrecursor>();
                foreach (double mz in spikes.Keys)
                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(mz, mixedPrecursor.MZ, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                        foreach (Sample sample in spikes[mz].Keys)
                            Isomers.Add(sample, spikes[mz][sample]);

                Dictionary<double, int> dicOfCommonFragments = GetCommonFragmentMz(Isomers.Values, nbProductsToKeep);

                Dictionary<Sample, List<ProductMatch>> matches = new Dictionary<Sample, List<ProductMatch>>();
                foreach (Sample sample in Isomers.Keys)
                    matches.Add(sample, Isomers[sample].Psms.GetCombinedSpectrum(dbOptions, dicOfCommonFragments));

                Dictionary<Sample, ElutionCurve> curves = new Dictionary<Sample,ElutionCurve>();
                foreach (Query query in mixedPrecursor.Queries)
                {
                    double timeInMilliSeconds = query.spectrum.RetentionTimeInMin * 60.0 * 1000.0;
                    double overFlow = 0;
                    double underFlow = 0;
                    double percentError = 0;
                    Dictionary<Sample, double> finalRatios = LaunchMaxFlowFromSpectrum(matches, 1000, query.spectrum.Peaks, dbOptions.productMassTolerance,
                                                            0, query.spectrum.PrecursorIntensityPerMilliSecond * query.spectrum.InjectionTime, ref overFlow, ref underFlow, ref percentError, dbOptions.ConSole);

                    if (percentError < 0.5)
                    {
                        //UsedPrecursorArea += localArea;
                        Dictionary<Sample, double> avgQuantifiedRatios = new Dictionary<Sample, double>();

                        foreach (Sample sRatio in finalRatios.Keys)
                        {
                            if (!curves.ContainsKey(sRatio))
                                curves.Add(sRatio, new ElutionCurve());

                            curves[sRatio].AddPoint(timeInMilliSeconds, finalRatios[sRatio] * query.spectrum.PrecursorIntensityPerMilliSecond);
                        }
                    }
                    else
                        Console.WriteLine("Ignored Spectrum : " + percentError);
                }//End of foreach query
                foreach(Sample sRatio in curves.Keys)
                {                    
                    curves[sRatio].Compute();
                    if (!values.ContainsKey(sRatio))
                        values.Add(sRatio, curves[sRatio].Area);
                    else
                        values[sRatio] += curves[sRatio].Area;
                }
            }//End of for each nbProduct            

            Dictionary<Sample, double> averagedValues = new Dictionary<Sample, double>();
            foreach(Sample sRatio in values.Keys)
                averagedValues.Add(sRatio, values[sRatio] / (double)(1 + maxProducts - minProducts));

            return averagedValues;
        }

        private static double ComputeMaxFlow(List<List<ProductMatch>> spikedMatches,
                                    List<MsMsPeak> mixedSpectrum, MassTolerance tolerance,
                                ref List<List<double>> optimalSolutions,
                                ref double percentError,
                                ref List<long> average, IConSol ConSole)
        {
            //Create dictionnary of usefull peaks
            Dictionary<float, double> mixedFragDic = new Dictionary<float, double>();
            foreach (List<ProductMatch> fragmentRatio in spikedMatches)
            {
                foreach (ProductMatch match in fragmentRatio)
                {
                    if (!mixedFragDic.ContainsKey((float)match.theoMz))
                    {
                        float closest = -1;
                        foreach (float key in mixedFragDic.Keys)
                            if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(match.theoMz, key, tolerance.Units)) <= tolerance.Value)
                                closest = key;
                        if (closest > 0)
                        {
                            ConSole.WriteLine("Potential problem with selected fragment masses ");
                            match.theoMz = closest;
                        }
                        else
                            mixedFragDic.Add((float)match.theoMz, 0);
                    }
                }
            }

            //Fill dictionnary with Intensities
            List<float> keys = new List<float>(mixedFragDic.Keys);
            foreach (MsMsPeak peak in mixedSpectrum)                
                foreach (float mz in keys)
                {
                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(peak.MZ, mz, tolerance.Units)) <= tolerance.Value)
                        mixedFragDic[mz] += peak.Intensity;
                }

            List<long> localFlows = new List<long>();
            foreach (List<ProductMatch> fragmentRatio in spikedMatches)
                localFlows.Add(FindLocalMaximumFlow(fragmentRatio, mixedFragDic));

            Dictionary<float, double> virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
            double overError = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);
            double underError = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
            double[] bestIndexes = new double[spikedMatches.Count];

            int iterSize = 1;
            double bestOverallError = double.MaxValue;
            List<long> bestLocalFlows = new List<long>();
            Random rnd = new Random();
            while (overError > 1 && iterSize < 10000)//anything less than 1 is an acceptable solution
            {
                for (int index = 0; index < bestIndexes.Length; index++)
                    bestIndexes[index] = -1;
                
                for (int i = 0; i < spikedMatches.Count; i++)
                {
                    if (localFlows[i] > 0)
                    {
                        localFlows[i] -= iterSize;
                        virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
                        double tmpErrorMinus = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
                        double tmpErrorPlus = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);

                        double tmpFlowRate = Math.Abs(overError - tmpErrorPlus);
                        double underDiff = Math.Abs(underError - tmpErrorMinus);
                        if (underDiff >= 1)
                            tmpFlowRate /= underDiff;
                        bestIndexes[i] = tmpFlowRate;
                        //if (double.IsNaN(tmpFlowRate) || double.IsInfinity(tmpFlowRate))
                        //    ConSole.WriteLine("schnit");
                        /*
                        if (tmpFlowRate == worstFlowRate)
                            Console.WriteLine("testse");
                        if(tmpFlowRate > worstFlowRate)
                        //if (tmpErrorPlus < overError && (tmpErrorMinus < smallestUnderError
                        //    || (tmpErrorMinus == smallestUnderError && tmpErrorPlus < smallestOverError)))
                        {
                            worstFlowRate = tmpFlowRate;
                            //smallestOverError = tmpErrorPlus;
                            //smallestUnderError = tmpErrorMinus;
                            bestIndexes[i] = tmpFlowRate;
                        }//*/
                        localFlows[i] += iterSize;
                    }
                }

                //Pick pseudo randomly best index
                double worstFlowRate = 0.0;
                for (int index = 0; index < bestIndexes.Length; index++)
                    if (bestIndexes[index] > worstFlowRate)
                    {
                        worstFlowRate = bestIndexes[index];
                    }

                if (worstFlowRate > 0)
                {
                    int nbMatching = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if(bestIndexes[index] >= worstFlowRate)
                            nbMatching++;

                    int iterChoice = rnd.Next(0, nbMatching - 1);
                    int iterNb = 0;
                    for (int index = 0; index < bestIndexes.Length; index++)
                        if (bestIndexes[index] >= worstFlowRate)
                        {
                            if (iterChoice == iterNb)
                                localFlows[index] -= iterSize;
                            iterNb++;
                        }
                    iterSize = 1;
                }
                else
                    iterSize++;

                virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);
                overError = MaxFlowHelper.ComputeOverflow(virtualSpectrum, mixedFragDic);
                underError = MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
                if (overError + underError < bestOverallError)
                {
                    bestLocalFlows = new List<long>(localFlows);
                    bestOverallError = overError + underError;
                }
            }
            optimalSolutions.Clear();

            List<double> newList = new List<double>();
            foreach (long localFlow in localFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);

            newList = new List<double>();
            foreach (long localFlow in bestLocalFlows)
                newList.Add(localFlow);
            optimalSolutions.Add(newList);
                                    
            //Compute average
            if (bestOverallError < double.MaxValue)
            {
                average.Clear();
                for (int i = 0; i < optimalSolutions[0].Count; i++)
                {
                    double sum = 0.0;
                    foreach (List<double> solution in optimalSolutions)
                        sum += solution[i];
                    double avg = sum / (double)optimalSolutions.Count;
                    average.Add((long)avg);
                }
            }

            //Compute expected error in percentage
                double sumOfIntensities = 0;
            foreach (double val in mixedFragDic.Values)
                sumOfIntensities += val;
            percentError = underError / sumOfIntensities;

            virtualSpectrum = BuildVirtualSpectrum(spikedMatches, localFlows, mixedFragDic);            

            return MaxFlowHelper.ComputeUnderflow(virtualSpectrum, mixedFragDic);
        }//*/

        private static Dictionary<float, double> BuildVirtualSpectrum(List<List<ProductMatch>> fragRatios, List<long> ratios, Dictionary<float, double> fragments)
        {
            Dictionary<float, double> virtualFrag = new Dictionary<float, double>(fragments.Count);
            for(int i = 0; i < ratios.Count; i++)
            {
                foreach (ProductMatch match in fragRatios[i])
                {
                    if (!virtualFrag.ContainsKey((float)match.theoMz))
                        virtualFrag.Add((float)match.theoMz, match.obsIntensity * ratios[i]);
                    else
                        virtualFrag[(float)match.theoMz] += match.obsIntensity * ratios[i];
                }
            }
            return virtualFrag;
        }

        private static long FindLocalMaximumFlow(List<ProductMatch> spikedMatches, Dictionary<float, double> fragments)
        {
            Dictionary<float, double> virtualFrag = new Dictionary<float, double>(fragments.Count);
            
            bool letsGo = false;
            double maxCumul = 0;
            foreach (ProductMatch match in spikedMatches)
            {
                if (!virtualFrag.ContainsKey((float)match.theoMz))
                    virtualFrag.Add((float)match.theoMz, match.obsIntensity);
                else
                    virtualFrag[(float)match.theoMz] += match.obsIntensity;
                if(virtualFrag[(float)match.theoMz] > 0)
                    letsGo = true;
                maxCumul += match.obsIntensity;
            }
            maxCumul *= 10000000;
            long nbCumul = 0;
            if (letsGo)
            {
                nbCumul = 1;
                while (nbCumul < maxCumul)
                {
                    foreach (float key in fragments.Keys)
                    {
                        if (virtualFrag.ContainsKey(key))
                        {
                            if (virtualFrag[key] * nbCumul > fragments[key] + 1)
                                return nbCumul;
                        }                 
                    }
                    nbCumul++;
                }
            }
            return nbCumul;
        }

        private static Dictionary<Sample, double> LaunchMaxFlowFromSpectrum(Dictionary<Sample, List<ProductMatch>> ratiosToFit,  
                                            int precision, GraphML_List<MsMsPeak> capacity, MassTolerance tolerance, 
                                            int returnType,//0 for max flow, 1 for best flow, 2 for average
                                            double PrecursorIntensityInCTrap,
                                            ref double overFlow, ref double underFlow, ref double errorInPercent, IConSol ConSole)
        {
            List<List<double>> solutions = new List<List<double>>();
            List<long> average = new List<long>();
            List<MsMsPeak> expandedCapacity = new List<MsMsPeak>();
            double sumOfProducts = 0;
            foreach (MsMsPeak peak in capacity)
            {
                double intensity = peak.Intensity;// +peak.Intensity * boostRatioForSpectrum;
                expandedCapacity.Add(new MsMsPeak(peak.MZ, intensity * precision, peak.Charge));
                sumOfProducts += peak.Intensity;
            }
            List<List<ProductMatch>> tmpRatiosToFit = new List<List<ProductMatch>>();
            foreach (List<ProductMatch> list in ratiosToFit.Values)
            {
                List<ProductMatch> pms = new List<ProductMatch>();
                foreach (ProductMatch pm in list)
                {
                    ProductMatch newPm = new ProductMatch(pm);
                    //newPm.obsIntensity *= sumOfProducts;
                    newPm.obsIntensity = newPm.normalizedIntensity * PrecursorIntensityInCTrap;
                    pms.Add(newPm);
                }
                tmpRatiosToFit.Add(pms);
            }

            double error = ComputeMaxFlow(tmpRatiosToFit, expandedCapacity, tolerance, ref solutions, ref errorInPercent, ref average, ConSole);

            double sumOfIntensities = 0;
            foreach (MsMsPeak peak in expandedCapacity)
                sumOfIntensities += peak.Intensity;

            overFlow = 0;
            underFlow = error;

            List<double> result = null;
            switch (returnType)
            {
                case 0:
                    result = GetResult(solutions[0], precision, underFlow, sumOfIntensities);
                    break;
                case 1:
                    result = GetResult(solutions[1], precision, underFlow, sumOfIntensities);
                    break;
                case 2:
                    List<double> tmpAverage = new List<double>();                    
                    foreach (double val in average)
                        tmpAverage.Add(val);
                    result = GetResult(tmpAverage, precision, underFlow, sumOfIntensities);
                    break;
            }
            Dictionary<Sample, double> resultPerSample = new Dictionary<Sample, double>();
            int i = 0;
            foreach (Sample key in ratiosToFit.Keys)
            {
                resultPerSample.Add(key, result[i]);
                i++;
            }
            return resultPerSample;
        }

        private static List<double> GetResult(List<double> solution, int precision, double underFlow, double sumOfIntensities)
        {
            List<double> rez = new List<double>();
            double sumVal = 0.0;
            foreach(double val in solution)
                sumVal += val;
            
            sumVal += (underFlow / sumOfIntensities) * precision;

            foreach (double val in solution)
                rez.Add(val / sumVal);

            return rez;
        }

        private static void BuildSinglePeptideVirtualSpectrum(Result precomputedResults, bool smoothPrecursor, 
                    int nbProductsToKeep, Dictionary<Sample, double> RatioNormalizer, ref Dictionary<double, Dictionary<Sample, List<ProductMatch>>> FinalSpikedProducts, 
                    ref Dictionary<double, Dictionary<Sample, ElutionCurve>> PrecursorCurves, ref Dictionary<double, Dictionary<Sample, double>> Normalizor, ref List<double> FragmentMz, int charge)
        {
            DBOptions dbOptions = precomputedResults.dbOptions;
            Samples Project = precomputedResults.samples;

            foreach (Sample sample in Project)
            {
                double peptideMass = Peptide.ComputeMonoisotopicMass(sample.nameColumn);

                double foundKey = 0.0;
                foreach (double key in FinalSpikedProducts.Keys)
                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(peptideMass, key, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                        foundKey = key;

                if (foundKey == 0.0)
                {
                    foundKey = peptideMass;
                    FinalSpikedProducts.Add(foundKey, new Dictionary<Sample,List<ProductMatch>>());
                    PrecursorCurves.Add(foundKey, new Dictionary<Sample, ElutionCurve>());
                    Normalizor.Add(foundKey, new Dictionary<Sample, double>());
                }
            }
            
            Dictionary<double, int> AllFragments = new Dictionary<double,int>();
            vsCSVWriter writer = new vsCSVWriter(dbOptions.OutputFolder + "FragmentsUsed_" + nbProductsToKeep + "Products.csv");
            writer.AddLine("Precursor Mass, Peptide Sequence,Fragment,Pos,Charge,Mz,Intensity,Normalized Intensity");
            foreach (double foundKey in FinalSpikedProducts.Keys)
            {
                Peptide peptide = null;
                Dictionary<Sample, List<ProductMatch>> SpikedProducts = new Dictionary<Sample,List<ProductMatch>>();
                //List<double> PeakAvgIntensities = new List<double>();
                //Get list of PSMs per sample
                double averagePrecursorArea = 0;
                //double avgMsMsProductIntensity = 0;
                //int nbAveragedSpectrum = 0;
                int nbAverageArea = 0;
                Dictionary<Sample, PeptideSpectrumMatches> listOfPSMs = new Dictionary<Sample, PeptideSpectrumMatches>();
                foreach (Sample sample in Project)
                {
                    //double avgPeakIntensity = 0;
                    PeptideSpectrumMatches psmList = new PeptideSpectrumMatches();
                    foreach (Query query in precomputedResults.queries)
                    {
                        if (query.precursor.Charge == charge && query.precursor.sample == sample)
                        {
                            foreach (PeptideSpectrumMatch psm in query.psms)
                            {
                                if (sample.nameColumn.CompareTo(psm.Peptide.Sequence) == 0)// && psm.MatchingProducts > 4)//.ProbabilityScore() > 0.020)// && psm.ProbabilityScore() > bestPsmScore)
                                {
                                    if (Math.Abs(Proteomics.Utilities.Numerics.CalculateMassError(psm.Peptide.MonoisotopicMass, foundKey, dbOptions.precursorMassTolerance.Units)) <= dbOptions.precursorMassTolerance.Value)
                                    {
                                        peptide = psm.Peptide;
                                        psmList.Add(psm);
                                        //avgPeakIntensity += query.spectrum.PrecursorIntensity;
                                        //avgMsMsProductIntensity += query.spectrum.TotalIntensity;
                                        //nbAveragedSpectrum++;
                                    }
                                }
                            }
                        }
                    }

                    if (psmList.Count > 0)
                    {
                        //avgPeakIntensity /= (double)psmList.Count;
                        psmList.Sort(PeptideSpectrumMatches.AscendingRetentionTime);
                        listOfPSMs.Add(sample, psmList);
                        List<ProductMatch> productList = psmList.GetCombinedSpectrum(precomputedResults.dbOptions);

                        productList.Sort(ProductMatch.AscendingWeightComparison);
                                                
                        nbAverageArea++;

                        SpikedProducts.Add(sample, productList);
                        PrecursorCurves[foundKey].Add(sample, ElutionCurve.Create(psmList));
                        averagePrecursorArea += PrecursorCurves[foundKey][sample].Area;
                        //PeakAvgIntensities.Add(avgPeakIntensity);
                    }
                    else
                    {
                        listOfPSMs.Add(sample, psmList);
                        PrecursorCurves[foundKey].Add(sample, ElutionCurve.Create(psmList));
                        //PeakAvgIntensities.Add(0);
                    }
                }//end of for each sample

                //avgMsMsProductIntensity /= nbAveragedSpectrum;
                averagePrecursorArea /= (double)nbAverageArea;

                Dictionary<double, int> DicOfFragmentsToKeep = new Dictionary<double, int>();
                //Get List of desired fragments, and keep masses
                foreach(Sample sample in SpikedProducts.Keys)
                {
                    for (int j = SpikedProducts[sample].Count - 1; j > 0 && j >= SpikedProducts[sample].Count - nbProductsToKeep; j--)
                    {
                        if (!DicOfFragmentsToKeep.ContainsKey(SpikedProducts[sample][j].theoMz))
                            DicOfFragmentsToKeep.Add(SpikedProducts[sample][j].theoMz, 0);

                        DicOfFragmentsToKeep[SpikedProducts[sample][j].theoMz]++;
                    }
                }

                foreach (Sample sample in SpikedProducts.Keys)
                {
                    List<ProductMatch> list = listOfPSMs[sample].GetCombinedSpectrum(precomputedResults.dbOptions, DicOfFragmentsToKeep);

                    foreach (ProductMatch match in list)
                        writer.AddLine(foundKey + "," + sample.nameColumn + "," + match.fragment + "," + match.fragmentPos + "," + match.charge + "," + match.theoMz + "," + match.obsIntensity + "," + match.normalizedIntensity);

                    FinalSpikedProducts[foundKey].Add(sample, list);
                }

                foreach(double key in DicOfFragmentsToKeep.Keys)
                    if(!AllFragments.ContainsKey(key))
                        AllFragments.Add(key, 1);
                    else
                        AllFragments[key]++;
                /*
                double avgPrecursors = 0;
                int nbInt = 0;
                foreach (double avg in PeakAvgIntensities)
                {
                    if (avg > 0)
                        nbInt++;
                    avgPrecursors += avg;
                }
                avgPrecursors /= (double)nbInt;
                /*
                for (int i = 0; i < Project.Count; i++)
                {
                    //Normalize each spectrum based on average precursor intensity
                    if (PeakAvgIntensities[i] > 0)
                    {
                        double y = (avgPrecursors - PeakAvgIntensities[i]) / PeakAvgIntensities[i];
                        for (int j = 0; j < FinalSpikedProducts[foundKey][i].Count; j++)
                            FinalSpikedProducts[foundKey][i][j].obsIntensity += y * FinalSpikedProducts[foundKey][i][j].obsIntensity;
                    }
                }//*/

                Dictionary<Sample, double> listOfsumOfProducts = new Dictionary<Sample, double>();
                double avgRatioSum = 0;
                int nbInt = 0;
                foreach (Sample sample in FinalSpikedProducts[foundKey].Keys)
                {
                    double sumOfProducts = 0;
                    for (int j = 0; j < FinalSpikedProducts[foundKey][sample].Count; j++)
                        sumOfProducts += FinalSpikedProducts[foundKey][sample][j].normalizedIntensity;
                    listOfsumOfProducts.Add(sample, sumOfProducts);

                    avgRatioSum += sumOfProducts; 

                    if (sumOfProducts > 0)
                        nbInt++;
                }
                avgRatioSum /= (double)nbInt;

                foreach (Sample sample in listOfsumOfProducts.Keys)
                {
                    if (listOfsumOfProducts[sample] > 0)
                    {
                        double lossOfPrecIntensity = 1;// Math.Log(averagePrecursorArea, 2) / Math.Log(PrecursorCurves[foundKey][sample].Area, 2);
                        double lossofFragmentIntensity = 1;// Math.Pow(avgRatioSum, 2) / Math.Pow(listOfsumOfProducts[sample], 2);

                        //Normalizor[foundKey].Add((spectrumWeight / listOfsumOfProducts[i]) * (averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        //Normalizor[foundKey].Add(lossOfPrecIntensity * lossofFragmentIntensity);//(averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        Normalizor[foundKey].Add(sample, lossofFragmentIntensity * lossOfPrecIntensity);//(averagePrecursorArea / PrecursorAreas[foundKey][i]) * RatioNormalizer[i]);
                        //for (int j = 0; j < FinalSpikedProducts[foundKey][i].Count; j++)
                        //    FinalSpikedProducts[foundKey][i][j].normalizedIntensity /= listOfsumOfProducts[i];
                    }
                    else
                        Normalizor[foundKey].Add(sample, 1.0);
                }
            }

            foreach(double key in AllFragments.Keys)
                FragmentMz.Add(key);

            writer.WriteToFile();
        }
    }
}
