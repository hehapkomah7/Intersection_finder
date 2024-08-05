using System;

namespace SegmentsIntersection
{
    public class Program
    {
        public static void Main()
        {
            Test();
        }

        static void Test()
        {
            string[] points =
            {
                "Enter separated by spaces coordinates X, Y, Z for the start of the first segment:",
                "Enter separated by spaces coordinates X, Y, Z for the end of the first segment:",
                "Enter separated by spaces coordinates X, Y, Z for the start of the second segment:",
                "Enter separated by spaces coordinates X, Y, Z for the end of the second segment:"
            };

            Vector3D[] vectors = new Vector3D[4];

            for (int i = 0; i < vectors.Length; i++)
            {
                Console.WriteLine(points[i]);
                string[] coordinates = Console.ReadLine().Split(' ');

                double x = double.Parse(coordinates[0]);
                double y = double.Parse(coordinates[1]);
                double z = double.Parse(coordinates[2]);

                vectors[i] = new Vector3D(x, y, z);
            }

            var segment1 = new Segment3D(vectors[0], vectors[1]);
            var segment2 = new Segment3D(vectors[2], vectors[3]);
            Intersector intersector = new Intersector();

            if(intersector.Checker(segment1, segment2))
            {
                Console.WriteLine("Segments are collinear or identical");
            }
            else
            {
                Vector3D point = intersector.Intersect(segment1, segment2);

                if (point == null)
                {
                    Console.WriteLine("Segments have no intersection");
                }
                else
                {
                    Console.WriteLine($"Intersection point: {point}");
                }
            }
        }
    }

    public class Vector3D
    {
        public double X;
        public double Y;
        public double Z;

        public Vector3D(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public override string ToString()
        {
            return $"({X}, {Y}, {Z})";
        }
    }

    public class Segment3D
    {
        public Vector3D Start;
        public Vector3D End;

        public Segment3D(Vector3D start, Vector3D end)
        {
            Start = start;
            End = end;
        }
    }

    public class Intersector
    {
        public bool Checker(Segment3D segment1, Segment3D segment2)
        {
            if (segment1.Start.Equals(segment2.Start) && segment1.End.Equals(segment2.End) || segment1.Start.Equals(segment2.End) && segment1.End.Equals(segment2.Start)) {
                return true;
            }

            Vector3D u = new Vector3D(segment1.End.X - segment1.Start.X, segment1.End.Y - segment1.Start.Y, segment1.End.Z - segment1.Start.Z);
            Vector3D v = new Vector3D(segment2.End.X - segment2.Start.X, segment2.End.Y - segment2.Start.Y, segment2.End.Z - segment2.Start.Z);
            Vector3D w = new Vector3D(segment1.Start.X - segment2.Start.X, segment1.Start.Y - segment2.Start.Y, segment1.Start.Z - segment2.Start.Z);

            if(CrossProduct(u, v).Value() < 1e-10 && CrossProduct(u, w).Value() < 1e-10)
            {
                return true;
            }

            return false;
        }

        public Vector3D Intersect(Segment3D segment1, Segment3D segment2)
        {
            Vector3D u = new Vector3D(segment1.End.X - segment1.Start.X, segment1.End.Y - segment1.Start.Y, segment1.End.Z - segment1.Start.Z);
            Vector3D v = new Vector3D(segment2.End.X - segment2.Start.X, segment2.End.Y - segment2.Start.Y, segment2.End.Z - segment2.Start.Z);
            Vector3D w = new Vector3D(segment1.Start.X - segment2.Start.X, segment1.Start.Y - segment2.Start.Y, segment1.Start.Z - segment2.Start.Z);

            double uu = DotProduct(u, u);
            double uv = DotProduct(u, v);
            double vv = DotProduct(v, v);
            double uw = DotProduct(u, w);
            double vw = DotProduct(v, w);

            //Console.WriteLine($"u: {u}, v: {v}, w: {w}");
            //Console.WriteLine($"uu: {uu}, uv: {uv}, vv: {vv}, uw: {uw}, vw: {vw}");

            double denominator = uu * vv - uv * uv;

            if (Math.Abs(denominator) < 1e-10) // Parallel or coincident segments
            {
                return null;
            }

            double s = (uv * vw - vv * uw) / denominator;
            double t = (uu * vw - uv * uw) / denominator;

            //Console.WriteLine($"s: {s}, t: {t}");

            if (s < 0 || s > 1 || t < 0 || t > 1)
            {
                return null; // Intersection is not within the segments
            }

            double X_sol = segment1.Start.X + s * u.X;
            double Y_sol = segment1.Start.Y + s * u.Y;
            double Z_sol = segment1.Start.Z + s * u.Z;

            return new Vector3D(X_sol, Y_sol, Z_sol);
        }

        private double DotProduct(Vector3D a, Vector3D b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        }

        private Vector3D CrossProduct(Vector3D a, Vector3D b)
        {
            return new Vector3D(
               a.Y * b.Z - a.Z * b.Y,
               a.Z * b.X - a.X * b.Z,
               a.X * b.Y - a.Y * b.X);
        }
    }

    public static class Vector3DExtensions
    {
        public static double Value(this Vector3D vector)
        {
            return Math.Sqrt(Math.Pow(vector.X, 2) + Math.Pow(vector.Y, 2) + Math.Pow(vector.Z, 2));
        }
    }
}

/*double[,] coeffitients = new double[3, 4];
coeffitients[0, 0] = 1 / p1;
coeffitients[0, 1] = -1 / p2;
coeffitients[0, 2] = 0;
coeffitients[0, 3] = segment1.Start.X / p1 - segment1.Start.Y / p2;
coeffitients[1, 0] = 0;
coeffitients[1, 1] = 1 / p2;
coeffitients[1, 2] = -1 / p3;
coeffitients[1, 3] = segment1.Start.Y / p2 - segment1.Start.Z / p3;
coeffitients[2, 0] = 0;
coeffitients[2, 1] = 1 / p5;
coeffitients[2, 2] = -1 / p6;
coeffitients[2, 3] = segment2.Start.Y / p5 - segment2.Start.Z / p6;
double v1 = coeffitients[0, 3];
double v2 = coeffitients[1, 3];
double v3 = coeffitients[2, 3];
double reversedDet = 1 / ((-1 / (p1 * p2 * p6)) + (1 / (p1 * p3 * p5)));
double X_sol = reversedDet * ((-1 / (p2 * p6) + 1 / (p3 * p5)) * v1 + (-1 / (p2 * p6)) * v2 + (1 / (p2 * p3)) * v3);
double Y_sol = reversedDet * ((-1 / (p1 * p6)) * v2 + (1 / (p1 * p3)) * v3);
double Z_sol = reversedDet * ((-1 / (p1 * p5)) * v2 + (1 / (p1 * p2)) * v3);
Vector3D result = new Vector3D(X_sol, Y_sol, Z_sol);

if (Math.Abs((X_sol / p4 - Y_sol / p5) - (segment2.Start.X / p4 - segment2.Start.Y / p5)) > 1e-9)
{
return null;
}

double Segment1_len = Math.Sqrt((Math.Pow(p1, 2) + Math.Pow(p2, 2) + Math.Pow(p3, 2)));
double Segment2_len = Math.Sqrt((Math.Pow(p4, 2) + Math.Pow(p5, 2) + Math.Pow(p6, 2)));
double Segment1_fp_len = Math.Sqrt(Math.Pow((X_sol - segment1.Start.X), 2) + Math.Pow((Y_sol - segment1.Start.Y), 2) + Math.Pow((Z_sol - segment1.Start.Z), 2));
double Segment1_sp_len = Math.Sqrt(Math.Pow((segment1.End.X - X_sol), 2) + Math.Pow((segment1.End.Y - Y_sol), 2) + Math.Pow((segment1.End.Z - Z_sol), 2));
double Segment2_fp_len = Math.Sqrt(Math.Pow((X_sol - segment2.Start.X), 2) + Math.Pow((Y_sol - segment2.Start.Y), 2) + Math.Pow((Z_sol - segment2.Start.Z), 2));
double Segment2_sp_len = Math.Sqrt(Math.Pow((segment2.End.X - X_sol), 2) + Math.Pow((segment2.End.Y - Y_sol), 2) + Math.Pow((segment2.End.Z - Z_sol), 2));

if (Math.Abs(Segment1_len - (Segment1_fp_len + Segment1_sp_len)) > 1e-9 && Math.Abs(Segment2_len - (Segment2_fp_len + Segment2_sp_len)) > 1e-9)
{
return null;
}

return result;*/

/*public double[] GaussMethod(double[,] matrix, double[] constants)
{
    int n = matrix.GetLength(0);
    double[,] augmentedMatrix = new double[n, n + 1];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            augmentedMatrix[i, j] = matrix[i, j];
        }
        augmentedMatrix[i, n] = constants[i];
    }

    // forward stroke
    for (int i = 0; i < n; i++)
    {
        // row with the maximum element in the column
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
        {
            if (Math.Abs(augmentedMatrix[k, i]) > Math.Abs(augmentedMatrix[maxRow, i]))
            {
                maxRow = k;
            }
        }

        // strings rearranging
        for (int k = i; k < n + 1; k++)
        {
            double temp = augmentedMatrix[i, k];
            augmentedMatrix[i, k] = augmentedMatrix[maxRow, k];
            augmentedMatrix[maxRow, k] = temp;
        }

        // triangular form
        for (int k = i + 1; k < n; k++)
        {
            double factor = augmentedMatrix[k, i] / augmentedMatrix[i, i];
            for (int j = i; j < n + 1; j++)
            {
                augmentedMatrix[k, j] -= factor * augmentedMatrix[i, j];
            }
        }
    }

    // reverse stroke
    double[] result = new double[n];
    for (int i = n - 1; i >= 0; i--)
    {
        result[i] = augmentedMatrix[i, n] / augmentedMatrix[i, i];
        for (int k = 0; k < i; k++)
        {
            augmentedMatrix[k, n] -= augmentedMatrix[k, i] * result[i];
        }
    }

    return result;
}*/

