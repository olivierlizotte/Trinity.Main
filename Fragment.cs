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
    /// <summary>
    /// Abstract class for typicall fragment behavior
    /// </summary>
    public abstract class FragmentClass : Modification
    {
        public abstract string Name { get; }
        public abstract IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul);
        public abstract bool IsReverse { get; }
        public abstract double Distribution { get; }
    }

    public class FragmentA : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "a"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return nTermCumul - 29.002741 + Constants.PROTON_MASS;//-CHO
        }
        override public bool IsReverse { get { return false; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentB : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "b"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            //yield return nTermCumul/* + 0;//*/ - 1.007276;//-H
            yield return nTermCumul + 0;//*/ - 1.007276;//-H
        }
        override public bool IsReverse { get { return false; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentC : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "c"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            //yield return nTermCumul/* + 17.02654915;//*/ +16.018724;//+NH2
            yield return nTermCumul + 17.02654915;//*/ +16.018724;//+NH2
        }
        override public bool IsReverse { get { return false; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentX : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "x"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            yield return cTermCumul + 43.9898346942;// +15.9949 + 12 + 18.0105646942 - 0.5;
            //yield return cTermCumul + 15.9949 + 12 + 18.0105646942 - 1.007276;
            //yield return cTermCumul + 26.98708959;//CO-H1.007276
        }
        override public bool IsReverse { get { return true; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentY : FragmentClass//Modification, IFragment 
    {
        public override string Name { get { return "y"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            //yield return cTermCumul/* + 18.0105646942;//*/ + 1.007276;//+H
            yield return cTermCumul + Constants.WATER_MONOISOTOPIC_MASS;//*/ + 1.007276;//+H
        }
        override public bool IsReverse { get { return true; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    public class FragmentZ : FragmentClass//Modification, IFragment
    {
        public override string Name { get { return "z"; } }
        override public IEnumerable<double> ComputeFragment(double cTermCumul, double nTermCumul)
        {
            //yield return cTermCumul/* + 1.991840552567;//*/ -16.018724;//-NH2
            yield return cTermCumul + 1.991840552567 - Constants.PROTON_MASS;//*/ -16.018724;//-NH2
        }
        override public bool IsReverse { get { return true; } }
        public override double Distribution
        {
            get { return 1; }
        }
    }

    /// <summary>
    /// Compute theoretical fragments based on a peptide sequence and the precursor charge
    /// </summary>
    public class Fragments : GraphML_List<FragmentClass>
    {
        public IEnumerable<ProductMatch> ComputeFragments(Peptide peptide, int precursorCharge, DBOptions dbOptions)
        {
            string sequence = peptide.BaseSequence;            
            ProductMatch match = new ProductMatch();

            double cumulativeNTerminalMass = 0;// 1.007276 + 15.9949;//TODO Add NTerminal modifications here            
            double cumulativeCTerminalMass = 0;// 1.007276;//TODO Add CTerminal modifications here            
            string cumulSeq = "";
            string cumulRevSeq = "";
            for (int r = 1; r <= sequence.Length; r++)
            {
                match.fragmentPos = r;
                //Peptide mods, computed during digestion
                //TODO Only search for modifications that were seen in the spectrum
                double tmp = 0;
                cumulSeq += sequence[r - 1];
                if (peptide.FixedModifications != null && peptide.FixedModifications.ContainsKey(r))
                    foreach (Modification mod in peptide.FixedModifications[r])
                        tmp += mod.MonoisotopicMassShift;
                if (peptide.VariableModifications != null && peptide.VariableModifications.ContainsKey(r))
                    tmp += peptide.VariableModifications[r].MonoisotopicMassShift;
                cumulativeNTerminalMass += AminoAcidMasses.GetMonoisotopicMass(sequence[r - 1]) + tmp;

                tmp = 0;
                cumulRevSeq += sequence[sequence.Length - r];
                if (peptide.FixedModifications != null && peptide.FixedModifications.ContainsKey(sequence.Length - r))
                    foreach (Modification mod in peptide.FixedModifications[sequence.Length - r])
                        tmp += mod.MonoisotopicMassShift;
                if (peptide.VariableModifications != null && peptide.VariableModifications.ContainsKey(sequence.Length - r))
                    tmp += peptide.VariableModifications[sequence.Length - r].MonoisotopicMassShift;
                cumulativeCTerminalMass += AminoAcidMasses.GetMonoisotopicMass(sequence[sequence.Length - r]) + tmp;

                for (int c = precursorCharge; c > 0; c--)
                {
                    match.charge = c;
                    foreach (FragmentClass fragment in this)
                    {
                        match.fragment = fragment.Name;
                        foreach (double product_mass in fragment.ComputeFragment(cumulativeCTerminalMass, cumulativeNTerminalMass))
                        {
                            match.weight = fragment.Distribution;//TODO Adjust this value by computing overall impact (times 10?)
                            //ModLess
                            match.theoMz = Numerics.MZFromMass(product_mass, c);
                            yield return match;

                            //Added mods
                            if (dbOptions.addFragmentMods)
                            {
                                foreach (Modification mod in FragmentDictionary.Fragments.Values)
                                {
                                    if ((fragment.IsReverse ? cumulRevSeq : cumulSeq).Contains(mod.AminoAcid))
                                    {
                                        match.fragment = mod.Description;
                                        match.weight = fragment.Distribution * mod.Probability;
                                        match.theoMz = Numerics.MZFromMass(product_mass + mod.MonoisotopicMassShift, c);
                                        yield return match;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (dbOptions.addFragmentLoss)
            {
                match.charge = 1;
                match.fragmentPos = -1;
                foreach (Modification mod in FragmentDictionary.AAFragments.Values)
                {
                    if (peptide.BaseSequence.Contains(mod.AminoAcid))
                    {
                        match.fragment = mod.Description;
                        match.theoMz = mod.MonoisotopicMassShift;
                        match.weight = mod.Probability;
                        yield return match;
                    }
                }
            }
        }
    }
}
