/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.Collections.Generic;

namespace Trinity.MaxQuant
{
    public static class ArrayUtil
    {
        private const int bootstrapBufferLen = 999;
        private const int maxBootstrapVectorLen = 99;
        public const int nBoots = 150;
        private static readonly int[,][] bootstrapBuffer = new int[bootstrapBufferLen, maxBootstrapVectorLen][];
        private static int count;
        public static Random random = new Random();

        static ArrayUtil()
        {
            for (int i = 0; i < bootstrapBufferLen; i++)
            {
                for (int j = 0; j < maxBootstrapVectorLen; j++)
                {
                    bootstrapBuffer[i, j] = new int[j];
                    for (int k = 0; k < j; k++)
                    {
                        bootstrapBuffer[i, j][k] = random.Next(j);
                    }
                }
            }
        }

        public static int[] GetBootstrapIndices(int n)
        {
            if (n < maxBootstrapVectorLen)
            {
                count = (count + 1) % bootstrapBufferLen;
                return bootstrapBuffer[count, n];
            }
            int[] result = new int[n];
            for (int i = 0; i < n; i++)
            {
                result[i] = random.Next(n);
            }
            return result;
        }

        /// <summary>
        /// Extracts the indexed elements from the given array.
        /// </summary>
        /// <typeparam name="T">Arbitrary type of the array elements.</typeparam>
        /// <param name="array">The input array.</param>
        /// <param name="indices">Indices of the elements to be extracted.</param>
        /// <returns>An array containing the elements of the input array indexed by the <code>indices</code> array.</returns>
        public static T[] SubArray<T>(T[] array, int[] indices)
        {
            T[] result = new T[indices.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = array[indices[i]];
            }
            return result;
        }

        /// <summary>
        /// Extracts the indexed elements from the given array.
        /// </summary>
        /// <typeparam name="T">Arbitrary type of the array elements.</typeparam>
        /// <param name="array">The input array.</param>
        /// <param name="indices">Indices of the elements to be extracted.</param>
        /// <returns>An array containing the elements of the input array indexed by the <code>indices</code> array.</returns>
        public static T[] SubArray<T>(T[] array, ushort[] indices)
        {
            T[] result = new T[indices.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = array[indices[i]];
            }
            return result;
        }

        /// <summary>
        /// Extracts the indexed elements from the given <code>List</code>.
        /// </summary>
        /// <typeparam name="T">Arbitrary type of the array elements.</typeparam>
        /// <param name="array">The input <code>List</code>.</param>
        /// <param name="indices">Indices of the elements to be extracted.</param>
        /// <returns>An array containing the elements of the input <code>List</code> indexed by the <code>indices</code> array.</returns>
        public static T[] SubArray<T>(List<T> array, int[] indices)
        {
            T[] result = new T[indices.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = array[indices[i]];
            }
            return result;
        }

        /// <summary>
        /// Extracts the first <code>len</code> elements from the input array. 
        /// </summary>
        /// <typeparam name="T">Arbitrary type of the array elements.</typeparam>
        /// <param name="array">The input array.</param>
        /// <param name="len">Length of the output array.</param>
        /// <returns>The first <code>len</code> elements of the input array.</returns>
        public static T[] SubArray<T>(T[] array, int len)
        {
            if (array.Length <= len)
            {
                return (T[])array.Clone();
            }
            T[] result = new T[len];
            for (int i = 0; i < len; i++)
            {
                result[i] = array[i];
            }
            return result;
        }

        public static int FloorIndex<T>(T[] array, T value) where T : IComparable<T>
        {
            int n = array.Length;
            if (n == 0)
            {
                return -1;
            }
            if (value.CompareTo(array[n - 1]) > 0)
            {
                return n - 1;
            }
            if (value.CompareTo(array[0]) < 0)
            {
                return -1;
            }
            int a = Array.BinarySearch(array, value);
            if (a >= 0)
            {
                return a;
            }
            return -2 - a;
        }

        public static int CeilIndex<T>(T[] array, T value) where T : IComparable<T>
        {
            int n = array.Length;
            if (n == 0)
            {
                return -1;
            }
            if (value.CompareTo(array[n - 1]) > 0)
            {
                return -1;
            }
            if (value.CompareTo(array[0]) < 0)
            {
                return 0;
            }
            int a = Array.BinarySearch(array, value);
            if (a >= 0)
            {
                return a;
            }
            return -1 - a;
        }

        public static int MaxInd(double[] x)
        {
            int n = x.Length;
            double max = double.MinValue;
            int ind = -1;
            for (int i = 0; i < n; i++)
            {
                double val = x[i];
                if (val > max)
                {
                    max = val;
                    ind = i;
                }
            }
            return ind;
        }

        public static double Median(double[] x)
        {
            int n = x.Length;
            if (n == 0)
            {
                return double.NaN;
            }
            int[] o = Order(x);
            if (n % 2 == 1)
            {
                return x[o[n / 2]];
            }
            return 0.5 * (x[o[n / 2 - 1]] + x[o[n / 2]]);
        }

        public static int[] Complement(int[] w, int n)
        {
            HashSet<int> dummy = new HashSet<int>(w);
            //dummy.AddAll(w);
            List<int> result = new List<int>();
            for (int i = 0; i < n; i++)
            {
                if (!dummy.Contains(i))
                {
                    result.Add(i);
                }
            }
            return result.ToArray();
        }

        public static int IndexOf<T>(T[] p, T q)
        {
            for (int i = 0; i < p.Length; i++)
            {
                if (p[i].Equals(q))
                {
                    return i;
                }
            }
            return -1;
        }

        public static double Sum(double[] x)
        {
            int n = x.Length;
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += x[i];
            }
            return sum;
        }

        public static float Sum(float[] x)
        {
            int n = x.Length;
            float sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += x[i];
            }
            return sum;
        }

        public static int Sum(int[] x)
        {
            int n = x.Length;
            int sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += x[i];
            }
            return sum;
        }

        public static double Variance(double[] x)
        {
            int n = x.Length;
            double mean = Mean(x);
            double var = 0;
            for (int i = 0; i < n; i++)
            {
                double w = x[i] - mean;
                var += w * w;
            }
            var /= (n - 1);
            return var;
        }

        public static double StandardDeviation(double[] x)
        {
            return Math.Sqrt(Variance(x));
        }


        public static T[] SubArray<T>(T[] array, int startIndex, int stopIndex)
        {
            int len = stopIndex - startIndex;
            T[] result = new T[len];
            for (int i = 0; i < len; i++)
            {
                result[i] = array[startIndex + i];
            }
            return result;
        }

        public static List<T> SubList<T>(List<T> list, int[] indices)
        {
            List<T> result = new List<T>();
            foreach (int index in indices)
            {
                result.Add(list[index]);
            }
            return result;
        }

        public static int[] ConsecutiveInts(int from, int to)
        {
            int len = to - from;
            int[] result = new int[len];
            for (int i = 0; i < len; i++)
            {
                result[i] = from + i;
            }
            return result;
        }

        public static T[] Concat<T>(T[] a, T[] b)
        {
            T[] result = new T[a.Length + b.Length];
            Array.Copy(a, 0, result, 0, a.Length);
            Array.Copy(b, 0, result, a.Length, b.Length);
            return result;
        }

        public static T[] Concat<T>(T[] a, T b)
        {
            T[] result = new T[a.Length + 1];
            Array.Copy(a, 0, result, 0, a.Length);
            result[a.Length] = b;
            return result;
        }

        public static T[] Concat<T>(T[][] x)
        {
            int len = 0;
            for (int i = 0; i < x.Length; i++)
            {
                len += x[i].Length;
            }
            T[] result = new T[len];
            int c = 0;
            for (int i = 0; i < x.Length; i++)
            {
                for (int j = 0; j < x[i].Length; j++)
                {
                    result[c++] = x[i][j];
                }
            }
            return result;
        }

        public static T[] UniqueValues<T>(T[] array) where T : IComparable<T>
        {
            if (array.Length < 2)
            {
                return array;
            }
            T[] sorted = (T[])array.Clone();
            Array.Sort(sorted);
            int counter = 1;
            T lastVal = sorted[0];
            for (int i = 1; i < sorted.Length; i++)
            {
                if (!lastVal.Equals(sorted[i]))
                {
                    lastVal = sorted[i];
                    sorted[counter++] = lastVal;
                }
            }
            return SubArray(sorted, counter);
        }

        public static double Min(double[] x)
        {
            int n = x.Length;
            if (n == 0)
            {
                return Double.NaN;
            }
            double min = Double.MaxValue;
            for (int i = 0; i < n; i++)
            {
                double val = x[i];
                if (val < min)
                {
                    min = val;
                }
            }
            return min;
        }

        public static int Min(int[] x)
        {
            int n = x.Length;
            if (n == 0)
            {
                return Int32.MaxValue;
            }
            int min = Int32.MaxValue;
            for (int i = 0; i < n; i++)
            {
                int val = x[i];
                if (val < min)
                {
                    min = val;
                }
            }
            return min;
        }

        public static int MinInd(double[] x)
        {
            int n = x.Length;
            double min = Double.MaxValue;
            int ind = -1;
            for (int i = 0; i < n; i++)
            {
                double val = x[i];
                if (val < min)
                {
                    min = val;
                    ind = i;
                }
            }
            return ind;
        }

        public static double Max(double[] x)
        {
            int n = x.Length;
            if (n == 0)
            {
                return Double.NaN;
            }
            double max = -Double.MaxValue;
            for (int i = 0; i < n; i++)
            {
                double val = x[i];
                if (val > max)
                {
                    max = val;
                }
            }
            return max;
        }

        public static float Max(float[] x)
        {
            int n = x.Length;
            if (n == 0)
            {
                return float.NaN;
            }
            float max = -float.MaxValue;
            for (int i = 0; i < n; i++)
            {
                float val = x[i];
                if (val > max)
                {
                    max = val;
                }
            }
            return max;
        }

        public static int Max(int[] x)
        {
            int n = x.Length;
            if (n == 0)
            {
                return Int32.MinValue;
            }
            int max = Int32.MinValue;
            for (int i = 0; i < n; i++)
            {
                int val = x[i];
                if (val > max)
                {
                    max = val;
                }
            }
            return max;
        }

        public static int MaxInd(float[] x)
        {
            int n = x.Length;
            float max = float.MinValue;
            int ind = -1;
            for (int i = 0; i < n; i++)
            {
                float val = x[i];
                if (val > max)
                {
                    max = val;
                    ind = i;
                }
            }
            return ind;
        }

        public static double Correlation(double[] x, double[] y)
        {
            if (x.Length < 3)
            {
                return 0;
            }
            double mx = Mean(x);
            double my = Mean(y);
            double xx = 0;
            double yy = 0;
            double xy = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double wx = x[i] - mx;
                double wy = y[i] - my;
                xx += wx * wx;
                yy += wy * wy;
                xy += wx * wy;
            }
            double denom = xx * yy;
            if (denom > 0.0)
            {
                return xy / Math.Sqrt(denom);
            }
            return 0f;
        }

        public static double Cosine(double[] x, double[] y)
        {
            if (x.Length < 3)
            {
                return 0;
            }
            double xx = 0;
            double yy = 0;
            double xy = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double wx = x[i];
                double wy = y[i];
                xx += wx * wx;
                yy += wy * wy;
                xy += wx * wy;
            }
            double denom = xx * yy;
            if (denom > 0.0)
            {
                return xy / Math.Sqrt(denom);
            }
            return 0f;
        }

        public static double Mean(double[] x)
        {
            int n = x.Length;
            if (n == 0)
            {
                return Double.NaN;
            }
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += x[i];
            }
            return sum / n;
        }

        public static float[] SmoothMean(float[] m, int width)
        {
            float[] result = new float[m.Length];
            for (int i = 0; i < result.Length; i++)
            {
                int min = Math.Max(0, i - width);
                int max = Math.Min(result.Length - 1, i + width);
                result[i] = Average(m, min, max);
            }
            return result;
        }

        public static float[] SmoothMedian(float[] m, int width)
        {
            float[] result = new float[m.Length];
            for (int i = 0; i < result.Length; i++)
            {
                int min = Math.Max(0, i - width);
                int max = Math.Min(result.Length - 1, i + width);
                result[i] = Median(m, min, max);
            }
            return result;
        }

        public static double[] SmoothMean(double[] m, int width)
        {
            double[] result = new double[m.Length];
            for (int i = 0; i < result.Length; i++)
            {
                int min = Math.Max(0, i - width);
                int max = Math.Min(result.Length - 1, i + width);
                result[i] = Average(m, min, max);
            }
            return result;
        }

        public static double[] SmoothMedian(double[] m, int width)
        {
            double[] result = new double[m.Length];
            for (int i = 0; i < result.Length; i++)
            {
                int min = Math.Max(0, i - width);
                int max = Math.Min(result.Length - 1, i + width);
                result[i] = Median(m, min, max);
            }
            return result;
        }

        public static int[] CalcLocalMinPositions(double[] y)
        {
            List<int> result = new List<int>();
            for (int i = 2; i < y.Length - 2; i++)
            {
                double m2 = y[i - 2];
                double m1 = y[i - 1];
                double x = y[i];
                double p1 = y[i + 1];
                double p2 = y[i + 2];
                if (IsMin(x, m1, p1, m2, p2))
                {
                    result.Add(i);
                }
            }
            int[] minPos = result.ToArray();
            double[] minY = SubArray(y, minPos);
            int[] o = Order(minY);
            minPos = SubArray(minPos, o);
            return minPos;
        }

        public static int[] CalcLocalMinPositions(float[] y)
        {
            List<int> result = new List<int>();
            for (int i = 2; i < y.Length - 2; i++)
            {
                float m2 = y[i - 2];
                float m1 = y[i - 1];
                float x = y[i];
                float p1 = y[i + 1];
                float p2 = y[i + 2];
                if (IsMin(x, m1, p1, m2, p2))
                {
                    result.Add(i);
                }
            }
            int[] minPos = result.ToArray();
            float[] minY = SubArray(y, minPos);
            int[] o = Order(minY);
            minPos = SubArray(minPos, o);
            return minPos;
        }

        /// <summary>
        /// For the sake of simplicity do all sorting tasks in the project always and ever with this method.
        /// </summary>
        /// <typeparam name="T">The array type has to inherit IComparable in order to have a 
        /// criterion to sort on.</typeparam>
        /// <param name="x">The input data to be sorted.</param>
        /// <returns>An array of indices such that if x is accessed with those indices the values are in 
        /// ascending (or to be more precise, non-decending) order.</returns>
        public static int[] Order<T>(T[] x) where T : IComparable<T>
        {
            int[] order = ConsecutiveInts(0, x.Length);
            const int low = 0;
            int high = order.Length - 1;
            int[] dummy = new int[order.Length];
            Array.Copy(order, dummy, order.Length);
            SortImpl(x, order, dummy, low, high);
            return order;
        }

        public static double[] Rank<T>(T[] data) where T : IComparable<T>
        {
            return Rank(data, true);
        }

        /// <summary>
        /// Calculates the rank of the given data. The lowest rank value is 0.
        /// The input array type must inherit IComparable.
        /// </summary>
        public static double[] Rank<T>(T[] data, bool tieCorrection) where T : IComparable<T>
        {
            int n = data.Length;
            double[] rank = new double[n];
            int[] index = Order(data);
            for (int j = 0; j < n; j++)
            {
                rank[index[j]] = j;
            }
            /* Fix for equal ranks */
            if (tieCorrection)
            {
                int i = 0;
                while (i < n)
                {
                    T value = data[index[i]];
                    int j = i + 1;
                    while (j < n && data[index[j]].Equals(value))
                    {
                        j++;
                    }
                    int m = j - i;
                    double v1 = rank[index[i]] + (m - 1) / 2.0;
                    for (j = i; j < i + m; j++)
                    {
                        rank[index[j]] = v1;
                    }
                    i += m;
                }
            }
            return rank;
        }

        public static float Cosine(float[] x, float[] y)
        {
            if (x.Length < 3)
            {
                return 0;
            }
            double xx = 0;
            double yy = 0;
            double xy = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double wx = x[i];
                double wy = y[i];
                xx += wx * wx;
                yy += wy * wy;
                xy += wx * wy;
            }
            double denom = xx * yy;
            if (denom > 0.0)
            {
                return (float)(xy / Math.Sqrt(denom));
            }
            return 0f;
        }

        /// <summary>
        /// Private class that implements the sorting algorithm.
        /// </summary>
        private static void SortImpl<T>(T[] data, int[] orderDest, int[] orderSrc, int low, int high)
            where T : IComparable<T>
        {
            if (low >= high)
            {
                return;
            }
            int mid = low + ((high - low) >> 1);
            SortImpl(data, orderSrc, orderDest, low, mid);
            SortImpl(data, orderSrc, orderDest, mid + 1, high);
            if (data[orderSrc[mid]].CompareTo(data[orderSrc[mid + 1]]) <= 0)
            {
                Array.Copy(orderSrc, low, orderDest, low, high - low + 1);
                return;
            }
            if (data[orderSrc[low]].CompareTo(data[orderSrc[high]]) > 0)
            {
                int m = (high - low) % 2 == 0 ? mid : mid + 1;
                Array.Copy(orderSrc, low, orderDest, m, mid - low + 1);
                Array.Copy(orderSrc, mid + 1, orderDest, low, high - mid);
                return;
            }
            int tLow = low;
            int tHigh = mid + 1;
            for (int i = low; i <= high; i++)
            {
                if ((tLow <= mid) &&
                    ((tHigh > high) || (data[orderSrc[tLow]]).CompareTo(data[orderSrc[tHigh]]) <= 0))
                {
                    orderDest[i] = orderSrc[tLow++];
                }
                else
                {
                    orderDest[i] = orderSrc[tHigh++];
                }
            }
        }

        private static bool IsMin(double x, double m1, double p1, double m2, double p2)
        {
            if (x < m1 && x < p1)
            {
                return true;
            }
            if (x == m1 && x < m2 && x < p1)
            {
                return true;
            }
            if (x < m1 && x == p1 && x < p2)
            {
                return true;
            }
            if (x < m2 && x == m1 && x == p1 && x < p2)
            {
                return true;
            }
            return false;
        }

        private static float Median(float[] m, int min, int max)
        {
            int len = max - min + 1;
            if (len == 1)
            {
                return m[min];
            }
            if (len == 2)
            {
                return 0.5f * (m[min] + m[max]);
            }
            if (len == 3)
            {
                float m1 = m[min];
                float m2 = m[min + 1];
                float m3 = m[min + 2];
                if (m1 <= m2 && m2 <= m3)
                {
                    return m2;
                }
                if (m2 <= m3 && m3 <= m1)
                {
                    return m3;
                }
                if (m3 <= m1 && m1 <= m2)
                {
                    return m1;
                }
                if (m3 <= m2 && m2 <= m1)
                {
                    return m2;
                }
                if (m2 <= m1 && m1 <= m3)
                {
                    return m1;
                }
                if (m1 <= m3 && m3 <= m2)
                {
                    return m3;
                }
            }
            float[] x = new float[len];
            for (int i = 0; i < len; i++)
            {
                x[i] = m[min + i];
            }
            Array.Sort(x);
            if (len % 2 == 0)
            {
                int w = len / 2;
                return 0.5f * (x[w - 1] + x[w]);
            }
            else
            {
                int w = len / 2;
                return x[w];
            }
        }

        private static float Average(float[] m, int min, int max)
        {
            float sum = 0;
            for (int i = min; i <= max; i++)
            {
                sum += m[i];
            }
            return sum / (max - min + 1);
        }

        private static double Median(double[] m, int min, int max)
        {
            int len = max - min + 1;
            if (len == 1)
            {
                return m[min];
            }
            if (len == 2)
            {
                return 0.5f * (m[min] + m[max]);
            }
            if (len == 3)
            {
                double m1 = m[min];
                double m2 = m[min + 1];
                double m3 = m[min + 2];
                if (m1 <= m2 && m2 <= m3)
                {
                    return m2;
                }
                if (m2 <= m3 && m3 <= m1)
                {
                    return m3;
                }
                if (m3 <= m1 && m1 <= m2)
                {
                    return m1;
                }
                if (m3 <= m2 && m2 <= m1)
                {
                    return m2;
                }
                if (m2 <= m1 && m1 <= m3)
                {
                    return m1;
                }
                if (m1 <= m3 && m3 <= m2)
                {
                    return m3;
                }
            }
            double[] x = new double[len];
            for (int i = 0; i < len; i++)
            {
                x[i] = m[min + i];
            }
            Array.Sort(x);
            if (len % 2 == 0)
            {
                int w = len / 2;
                return 0.5f * (x[w - 1] + x[w]);
            }
            else
            {
                int w = len / 2;
                return x[w];
            }
        }

        private static double Average(double[] m, int min, int max)
        {
            double sum = 0;
            for (int i = min; i <= max; i++)
            {
                sum += m[i];
            }
            return sum / (max - min + 1);
        }

        public static int ClosestIndex(double[] array, double value)
        {
            int n = array.Length;
            if (n == 0)
            {
                return -1;
            }
            if (value.CompareTo(array[n - 1]) > 0)
            {
                return n - 1;
            }
            if (value.CompareTo(array[0]) < 0)
            {
                return 0;
            }
            int a = Array.BinarySearch(array, value);
            if (a >= 0)
            {
                return a;
            }
            if (array[-1 - a] - value < value - array[-2 - a])
            {
                return -1 - a;
            }
            return -2 - a;
        }

        public static int ClosestIndex(float[] array, float value)
        {
            int n = array.Length;
            if (n == 0)
            {
                return -1;
            }
            if (value.CompareTo(array[n - 1]) > 0)
            {
                return n - 1;
            }
            if (value.CompareTo(array[0]) < 0)
            {
                return 0;
            }
            int a = Array.BinarySearch(array, value);
            if (a >= 0)
            {
                return a;
            }
            if (array[-1 - a] - value < value - array[-2 - a])
            {
                return -1 - a;
            }
            return -2 - a;
        }
    }
}