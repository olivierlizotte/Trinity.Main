﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Trinity
{
    public class ProductMatch : GraphML_Node
    {
        public double theoMz;
        public double obsMz;
        public int charge;
        public double obsIntensity;
        public double mass_diff;
        public string fragment;
        public int fragmentPos;
        public double weight;
        public ProductMatch()
        {
        }

        public ProductMatch(ProductMatch match)
        {
            this.theoMz = match.theoMz;
            this.obsMz = match.obsMz;
            this.obsIntensity = match.obsIntensity;
            this.mass_diff = match.mass_diff;
            this.charge = match.charge;
            this.fragment = match.fragment;
            this.fragmentPos = match.fragmentPos;
            this.weight = match.weight;
        }
        public static int AscendingWeightComparison(ProductMatch left, ProductMatch right)
        {
            if (left.weight == right.weight)
                return left.obsIntensity.CompareTo(right.obsIntensity);
            return left.weight.CompareTo(right.weight);
        }
        public override string ToString()
        {
            return "[" + obsMz + ";" + obsIntensity + "]";
        }
    }
}
