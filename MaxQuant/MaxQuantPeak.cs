/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;
using System.IO;
using Proteomics.Utilities;
using Proteomics.Utilities.MaxQuant;

namespace Trinity.MaxQuant
{
    public class Peak
    {
        private double[] centerMz;
        private int isotopeClusterIndex = -1;
        private byte[] massRange;
        private float[] maxMz;
        private float[] minMz;
        private float[] origIntensity;
        private int[] scanIndex;
        private float[] smoothIntensity;

        private double totalCenterMz;
        private float totalCenterMzError;
        private float totalIntensity;
        private float totalMaxMz;
        private float totalMinMz;

        public Peak(double[] centerMz, float[] minMz, float[] maxMz, float[] smoothIntensity, float[] origIntensity,
                    int[] scanIndex, byte[] massRange)
        {
            this.centerMz = centerMz;
            this.minMz = minMz;
            this.maxMz = maxMz;
            this.smoothIntensity = smoothIntensity;
            this.origIntensity = origIntensity;
            this.scanIndex = scanIndex;
            this.massRange = massRange;
        }

        public Peak(BinaryReader reader)
        {
            totalCenterMz = reader.ReadDouble();
            totalCenterMzError = reader.ReadSingle();
            totalIntensity = reader.ReadSingle();
            totalMinMz = reader.ReadSingle();
            totalMaxMz = reader.ReadSingle();
            isotopeClusterIndex = reader.ReadInt32();
            int count = reader.ReadInt32();
            centerMz = new double[count];
            minMz = new float[count];
            maxMz = new float[count];
            smoothIntensity = new float[count];
            origIntensity = new float[count];
            scanIndex = new int[count];
            massRange = new byte[count];
            for (int i = 0; i < count; i++)
            {
                centerMz[i] = reader.ReadDouble();
                minMz[i] = reader.ReadSingle();
                maxMz[i] = reader.ReadSingle();
                smoothIntensity[i] = reader.ReadSingle();
                origIntensity[i] = reader.ReadSingle();
                scanIndex[i] = reader.ReadInt32();
                massRange[i] = reader.ReadByte();
            }
        }

        public double CenterMass
        {
            get { return totalCenterMz; }
        }

        public int Count
        {
            get { return centerMz.Length; }
        }

        public int LastScanIndex
        {
            get
            {
                if (Count == 0)
                {
                    return -1;
                }
                return scanIndex[Count - 1];
            }
        }

        public int FirstScanIndex
        {
            get
            {
                if (Count == 0)
                {
                    return -1;
                }
                return scanIndex[0];
            }
        }

        public float CenterMassError
        {
            get { return totalCenterMzError; }
        }

        public float Intensity
        {
            get { return totalIntensity; }
        }

        public int IsotopeClusterIndex
        {
            get { return isotopeClusterIndex; }
            set { isotopeClusterIndex = value; }
        }

        public double[] CenterMz
        {
            get { return centerMz; }
        }

        public float[] OrigIntensity
        {
            get { return origIntensity; }
        }

        public byte[] MassRange
        {
            get { return massRange; }
        }

        public float GetMinTime(IRawFile rawFile)
        {
            return (float)rawFile.GetMS1TimeSpan(Math.Max(0, FirstScanIndex - 1))[0];
        }

        public float GetMaxTime(IRawFile rawFile)
        {
            return (float)rawFile.GetMS1TimeSpan(LastScanIndex)[1];
        }

        public int GetScanIndex(int index)
        {
            return scanIndex[index];
        }

        public float GetSmoothIntensity(int index)
        {
            return smoothIntensity[index];
        }

        public double GetOriginalIntensity(int index)
        {
            return origIntensity[index];
        }

        public double GetCenterMass(int index)
        {
            return centerMz[index];
        }

        public double GetMinMass(int index)
        {
            return minMz[index];
        }

        public double GetMaxMass(int index)
        {
            return maxMz[index];
        }

        public void Write(BinaryWriter writer)
        {
            writer.Write(totalCenterMz);
            writer.Write(totalCenterMzError);
            writer.Write(totalIntensity);
            writer.Write(totalMinMz);
            writer.Write(totalMaxMz);
            writer.Write(isotopeClusterIndex);
            writer.Write(Count);
            for (int i = 0; i < Count; i++)
            {
                writer.Write(centerMz[i]);
                writer.Write(minMz[i]);
                writer.Write(maxMz[i]);
                writer.Write(smoothIntensity[i]);
                writer.Write(origIntensity[i]);
                writer.Write(scanIndex[i]);
                writer.Write(massRange[i]);
            }
        }

        public void CalcAverages(bool maxIntensity)
        {
            totalMinMz = float.MaxValue;
            totalMaxMz = -float.MaxValue;
            totalIntensity = 0;
            for (int i = 0; i < Count; i++)
            {
                if (minMz[i] < totalMinMz)
                {
                    totalMinMz = minMz[i];
                }
                if (maxMz[i] > totalMaxMz)
                {
                    totalMaxMz = maxMz[i];
                }
                if (maxIntensity)
                {
                    if (origIntensity[i] > totalIntensity)
                    {
                        totalIntensity = origIntensity[i];
                    }
                }
                else
                {
                    totalIntensity += origIntensity[i];
                }
            }
            double[] m = new double[ArrayUtil.nBoots];
            for (int i = 0; i < ArrayUtil.nBoots; i++)
            {
                m[i] = CalcAverageMz(this, ArrayUtil.GetBootstrapIndices(Count), null, null, 0.0);
            }
            totalCenterMz = 0;
            for (int i = 0; i < ArrayUtil.nBoots; i++)
            {
                totalCenterMz += m[i];
            }
            totalCenterMz /= ArrayUtil.nBoots;
            totalCenterMzError = 0;
            for (int i = 0; i < ArrayUtil.nBoots; i++)
            {
                totalCenterMzError += (float)((totalCenterMz - m[i]) * (totalCenterMz - m[i]));
            }
            totalCenterMzError /= ArrayUtil.nBoots;
            totalCenterMzError = (float)Math.Sqrt(totalCenterMzError);
        }

        public static double CalcAverageMz(Peak peak, int[] indx, double[,] mzCalibrationPar,
                                           double[,] intensityCalibrationPar, double absoluteCalibration)
        {
            if (double.IsNaN(absoluteCalibration))
            {
                absoluteCalibration = 0;
            }
            double norm1 = 0;
            double m = 0;
            if (mzCalibrationPar == null)
            {
                foreach (int index0 in indx)
                {
                    norm1 += peak.origIntensity[index0];
                    double mz = peak.centerMz[index0];
                    m += peak.origIntensity[index0] * mz;
                }
                m /= norm1;
                return m;
            }
            foreach (int index0 in indx)
            {
                norm1 += peak.origIntensity[index0];
                double mz = peak.centerMz[index0] * (1 + 1e-6 * NumericalRecipes.RelativeCorrection(peak.centerMz[index0], mzCalibrationPar,
                                                                                           intensityCalibrationPar,
                                                                                           peak.origIntensity[index0],
                                                                                           peak.massRange[index0]) +
                                                     1e-6 * absoluteCalibration);
                m += peak.origIntensity[index0] * mz;
            }
            m /= norm1;
            return m;
        }
        /*
        public double CalcAverageMz(double[,] mzCalibrationPar, double[,] intensityCalibrationPar, double absoluteCalibration)
        {
            int[] ind = ArrayUtil.ConsecutiveInts(0, Count);
            return CalcAverageMz(this, ind, mzCalibrationPar, intensityCalibrationPar, absoluteCalibration);
        }//*/

        public void RemoveDoublePeaks(IRawFile rawFile, bool maxIntensity)
        {
            CalcAverages(maxIntensity);
            int len = LastScanIndex - FirstScanIndex + 1;
            List<int>[] lists = new List<int>[len];
            for (int i = 0; i < len; i++)
            {
                lists[i] = new List<int>();
            }
            for (int i = 0; i < Count; i++)
            {
                lists[scanIndex[i] - FirstScanIndex].Add(i);
            }
            List<int> valids = new List<int>();
            foreach (List<int> list in lists)
            {
                if (list.Count == 1)
                {
                    valids.Add(list[0]);
                }
                else if (list.Count > 1)
                {
                    double[] dm = new double[list.Count];
                    for (int i = 0; i < list.Count; i++)
                    {
                        dm[i] = Math.Abs(centerMz[list[i]] - totalCenterMz);
                    }
                    int ind = ArrayUtil.Order(dm)[0];
                    valids.Add(list[ind]);
                }
            }
            int[] v = valids.ToArray();
            centerMz = ArrayUtil.SubArray(centerMz, v);
            minMz = ArrayUtil.SubArray(minMz, v);
            maxMz = ArrayUtil.SubArray(maxMz, v);
            smoothIntensity = ArrayUtil.SubArray(smoothIntensity, v);
            origIntensity = ArrayUtil.SubArray(origIntensity, v);
            scanIndex = ArrayUtil.SubArray(scanIndex, v);
            CalcAverages(maxIntensity);
        }

        public void Smooth(bool maxIntensity)
        {
            minMz = ArrayUtil.SmoothMean(minMz, 2);
            maxMz = ArrayUtil.SmoothMean(maxMz, 2);
            smoothIntensity = ArrayUtil.SmoothMedian(origIntensity, 1);
            smoothIntensity = ArrayUtil.SmoothMean(smoothIntensity, 2);
            CalcAverages(maxIntensity);
        }

        private Peak[] SplitIntoPair(double valleyFactor, bool maxIntensity)
        {
            float[] intensities = smoothIntensity;
            int[] minPos = ArrayUtil.CalcLocalMinPositions(intensities);
            foreach (int pos in minPos)
            {
                double leftMax = GetLeftMax(pos, intensities);
                double rightMax = GetRightMax(pos, intensities);
                double smallMax = Math.Min(leftMax, rightMax);
                if (smallMax / intensities[pos] > valleyFactor)
                {
                    return SplitAt(pos, maxIntensity);
                }
            }
            return null;
        }

        private Peak[] SplitAt(int pos, bool maxIntensity)
        {
            return new Peak[] { SubPeak(0, pos, maxIntensity), SubPeak(pos + 1, Count, maxIntensity) };
        }

        private Peak SubPeak(int start, int end, bool maxIntensity)
        {
            Peak result = new Peak(
                ArrayUtil.SubArray(centerMz, start, end),
                ArrayUtil.SubArray(minMz, start, end),
                ArrayUtil.SubArray(maxMz, start, end),
                ArrayUtil.SubArray(smoothIntensity, start, end),
                ArrayUtil.SubArray(origIntensity, start, end),
                ArrayUtil.SubArray(scanIndex, start, end),
                ArrayUtil.SubArray(massRange, start, end)
                );
            result.CalcAverages(maxIntensity);
            return result;
        }

        public void CorrectMasses(double[] dmx, bool maxIntensity)
        {
            for (int i = 0; i < centerMz.Length; i++)
            {
                int ind = scanIndex[i];
                double dm = dmx[ind];
                centerMz[i] -= centerMz[i] * dm * 1e-6;
                minMz[i] = (float)(minMz[i] - minMz[i] * dm * 1e-6);
                maxMz[i] = (float)(maxMz[i] - maxMz[i] * dm * 1e-6);
            }
            CalcAverages(maxIntensity);
        }

        private static float GetLeftMax(int pos, float[] y)
        {
            float max = -float.MaxValue;
            for (int i = 0; i < pos; i++)
            {
                if (y[i] > max)
                {
                    max = y[i];
                }
            }
            return max;
        }

        private static float GetRightMax(int pos, float[] y)
        {
            float max = -float.MaxValue;
            for (int i = pos + 1; i < y.Length; i++)
            {
                if (y[i] > max)
                {
                    max = y[i];
                }
            }
            return max;
        }

        public Peak[] Decompose(double valleyFactor, bool split, bool maxIntensity)
        {
            if (!split)
            {
                return new Peak[] { this };
            }
            List<Peak> splitCandidates = new List<Peak>();
            List<Peak> result = new List<Peak>();
            splitCandidates.Add(this);
            while (splitCandidates.Count > 0)
            {
                Peak maybeSplit = splitCandidates[0];
                splitCandidates.Remove(maybeSplit);
                Peak[] w = maybeSplit.SplitIntoPair(valleyFactor, maxIntensity);
                if (w == null)
                {
                    result.Add(maybeSplit);
                }
                else
                {
                    splitCandidates.Add(w[0]);
                    splitCandidates.Add(w[1]);
                    maybeSplit.Dispose();
                }
            }
            return result.ToArray();
        }

        public int[] GetScanIndices()
        {
            return scanIndex;
        }

        public float[] GetOriginalIntensities()
        {
            return origIntensity;
        }

        public void Dispose()
        {
            centerMz = null;
            minMz = null;
            maxMz = null;
            smoothIntensity = null;
            origIntensity = null;
            scanIndex = null;
            massRange = null;
        }

        public bool Contains(double m, double t, double minTime, double maxTime)
        {
            return (m <= totalMaxMz && m >= totalMinMz && t <= maxTime && t >= minTime);
        }
    }
}