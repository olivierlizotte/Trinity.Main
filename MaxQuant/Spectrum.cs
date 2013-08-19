/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;

namespace Trinity.MaxQuant
{
    public class Spectrum
    {
        private List<double> intensities;
        private List<double> masses;

        public Spectrum(pwiz.CLI.msdata.BinaryData masses, pwiz.CLI.msdata.BinaryData intensities)
        {            
            this.masses = new List<double>(masses);
            this.intensities = new List<double>(intensities);
        }
        public Spectrum(List<double> masses, List<double> intensities)
        {
            this.masses = masses;
            this.intensities = intensities;
        }
        public Spectrum(double[] masses, double[] intensities)
        {
            this.masses = new List<double>(masses);
            this.intensities = new List<double>(intensities);
        }

        public int SizeEstimate
        {
            get { return 12 * Count + 16; }
        }

        public int Count
        {
            get { return masses.Count; }
        }

        public double MinMass
        {
            get { return GetMass(0); }
        }

        public double MaxMass
        {
            get { return GetMass(Count - 1); }
        }

        public int GetCeilIndex(double mass)
        {
            if (Count == 0)
            {
                return -1;
            }
            if (mass <= MinMass)
            {
                return 0;
            }
            if (mass > MaxMass)
            {
                return -1;
            }
            
            int index = masses.BinarySearch(mass);
            if (index >= 0)
            {
                return index;
            }
            index = -1 - index;
            return index;
        }

        public int GetFloorIndex(double mass)
        {
            if (Count == 0)
            {
                return -1;
            }
            if (mass >= MaxMass)
            {
                return Count - 1;
            }
            if (mass < MinMass)
            {
                return -1;
            }
            int index = masses.BinarySearch(mass);
            if (index >= 0)
            {
                return index;
            }
            index = -2 - index;
            return index;
        }

        public int GetClosestIndex(double mass, bool outOfRangeIsInvalid)
        {
            if (mass <= MinMass)
            {
                return outOfRangeIsInvalid ? -1 : 0;
            }
            if (mass >= MaxMass)
            {
                return outOfRangeIsInvalid ? -1 : Count - 1;
            }
            int index = masses.BinarySearch(mass);
            if (index >= 0)
            {
                return index;
            }
            index = -2 - index;
            if (Math.Abs(GetMass(index) - mass) < Math.Abs(GetMass(index + 1) - mass))
            {
                return index;
            }
            return index + 1;
        }

        public double GetMass(int index)
        {
            return masses[index];
        }

        public double GetIntensity(int index)
        {
            return intensities[index];
        }

        public double GetIntensityFromMass(double mass)
        {
            int ind = GetClosestIndex(mass, true);
            if (ind == -1)
            {
                return 0;
            }
            return GetIntensity(ind);
        }

        public bool IsMax(double x, double m1, double p1, double m2, double p2)
        {
            if (x > m1 && x > p1)
            {
                return true;
            }
            if (x > m2 && x == m1 && x > p1)
            {
                return true;
            }
            if (x > m1 && x == p1 && x > p2)
            {
                return true;
            }
            return false;
        }

        public int CalcMinPeakIndex(int ind)
        {
            while (ind > 0 && intensities[ind - 1] != 0 && intensities[ind - 1] < intensities[ind])
            {
                ind--;
            }
            return ind;
        }

        public int CalcMaxPeakIndex(int ind)
        {
            while (ind < Count - 1 && intensities[ind + 1] != 0 && intensities[ind + 1] < intensities[ind])
            {
                ind++;
            }
            return ind;
        }

        public void Dispose()
        {
            masses = null;
            intensities = null;
        }

        public double[] GetMasses()
        {
            double[] result = new double[Count];
            for (int i = 0; i < Count; i++)
            {
                result[i] = GetMass(i);
            }
            return result;
        }

        public double[] GetIntensities()
        {
            double[] result = new double[Count];
            for (int i = 0; i < Count; i++)
            {
                result[i] = GetIntensity(i);
            }
            return result;
        }

        public Spectrum TopX(int topX, double window)
        {
            int[] index = TopXIndices(topX, window);
            return new Spectrum(ArrayUtil.SubArray(masses, index), ArrayUtil.SubArray(intensities, index));
        }

        private int[] TopXIndices(int topx, double window)
        {
            int[] top = new int[topx];
            for (int n = 0; n < topx; n++)
            {
                top[n] = -1;
            }
            List<int> index = new List<int>();
            int leftSide = 0;
            for (int i = 0; i < Count; i++)
            {
                if (GetMass(i) - GetMass(leftSide) <= window)
                {
                    for (int j = 0; j < topx; j++)
                    {
                        if (top[j] == -1 || GetIntensity(i) >= GetIntensity(top[j]))
                        {
                            for (int k = topx - 1; k > j; k--)
                            {
                                top[k] = top[k - 1];
                            }
                            top[j] = i;
                            break;
                        }
                    }
                }
                else
                {
                    Array.Sort(top);
                    for (int id = 0; id < topx; id++)
                    {
                        if (top[id] >= 0)
                        {
                            index.Add(top[id]);
                        }
                    }
                    for (int n = 0; n < topx; n++)
                    {
                        top[n] = -1;
                    }
                    leftSide = i;
                    i--;
                }
                if (i + 1 == Count)
                {
                    Array.Sort(top);
                    for (int id = 0; id < topx; id++)
                    {
                        if (top[id] >= 0)
                        {
                            index.Add(top[id]);
                        }
                    }
                    for (int n = 0; n < topx; n++)
                    {
                        top[n] = -1;
                    }
                    leftSide = i;
                }
            }
            return index.ToArray();
        }
    }
}