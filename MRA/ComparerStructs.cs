using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MRA
{
    class DistinctDoubleArrayComparer : IEqualityComparer<double[]>
    {
        public bool Equals(double[] x, double[] y)
        {
            if (x.Length != y.Length) { return false; }
            else if (x.Length != 3 || y.Length != 3) { return false; }

            return x[0] == y[0] && x[1] == y[1] && x[2] == y[2];
        }

        public int GetHashCode(double[] obj)
        {
            return -1;
        }
    }

    public struct DataTempLoop
    {
        public double[] sCurrent;
        public double[] sAfter;
        public double[] lCurrent;
        public double[] lAfter;
        public double[,] uCurrent;
        public double[,] uAfter;
        public double[,] vCurrent;
        public double[,] vAfter;
        public double[,] aCurrent;
        public double[,] aAfter;
        public double[] norm_v;
        public double[] norm_a;
    }

    public struct FinalResults
    {
        public List<int> lasche1;
        public int numLoose;
        public double meanLoose;
        public double minLoose;
        public int numOver;
        public double meanOver;
        public double maxOver;
        public double perDiff;
        public int n_ropes;
        public List<double> maxLengths;
        public List<double> lFinal;
        public List<int[]> link;
        public double[,] uFinal;
    }

}
