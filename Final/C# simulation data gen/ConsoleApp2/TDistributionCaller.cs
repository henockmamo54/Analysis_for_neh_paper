using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp2
{

    public class TDistributionCaller
    {

        double[] w = new double[] { 0.0055657196642445571,
                                    0.012915947284065419,0.020181515297735382,0.027298621498568734,
                                    0.034213810770299537,0.040875750923643261,0.047235083490265582,
                                    0.053244713977759692,0.058860144245324798,0.064039797355015485,
                                    0.068745323835736408,0.072941885005653087,0.076598410645870640,
                                    0.079687828912071670,0.082187266704339706,0.084078218979661945,
                                    0.085346685739338721,0.085983275670394821 };

        double[] y = new double[]  {0.0021695375159141994,
                                    0.011413521097787704,0.027972308950302116,0.051727015600492421,
                                    0.082502225484340941, 0.12007019910960293,0.16415283300752470,
                                    0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
                                    0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
                                    0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
                                    0.87126389619061517, 0.95698180152629142 };
        double DBL_MIN = 2.2250738585072014e-308;


        double gammln(double xx)
        {
            int j;

            double x, tmp, y, ser;

            double[] cof = new double[]{ 57.1562356658629235,-59.5979603554754912,
                14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
                .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
                -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
                .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5 };

            if (xx <= 0)
                throw new Exception("bad arg in gammln");

            y = x = xx;

            tmp = x + 5.24218750000000000;

            tmp = (x + 0.5) * Math.Log(tmp) - tmp;

            ser = 0.999999999999997092;

            for (j = 0; j < 14; j++)
            {
                ser += cof[j] / ++y;
            }

            return tmp + Math.Log(2.5066282746310005 * ser / x);
        }


        double betacf(double a, double b, double x)
        {

            int m, m2;

            double aa, c, d, del, h, qab, qam, qap;

            double EPS = 1E-8; // 1.e - 8;

            double FPMIN = DBL_MIN / EPS;

            qab = a + b;

            qap = a + 1.0;

            qam = a - 1.0;

            c = 1.0;

            d = 1.0 - qab * x / qap;

            if (Math.Abs(d) < FPMIN) d = FPMIN;

            d = 1.0 / d;

            h = d;

            for (m = 1; m < 10000; m++)
            {
                m2 = 2 * m;
                aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = 1.0 + aa * d;
                if (Math.Abs(d) < FPMIN)
                    d = FPMIN;

                c = 1.0 + aa / c;

                if (Math.Abs(c) < FPMIN)
                    c = FPMIN;

                d = 1.0 / d;

                h *= d * c;

                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));

                d = 1.0 + aa * d;

                if (Math.Abs(d) < FPMIN)
                    d = FPMIN;

                c = 1.0 + aa / c;

                if (Math.Abs(c) < FPMIN)
                    c = FPMIN;

                d = 1.0 / d;

                del = d * c;

                h *= del;

                if (Math.Abs(del - 1.0) <= EPS)
                    break;
            }
            return h;

        }

        double betai(double a, double b, double x)
        {
            double bt;

            int SWITCH = 3000;

            if (a <= 0.0 || b <= 0.0)
            {
                throw new Exception("Bad a or b in routine betai");
            }

            if (x < 0.0 || x > 1.0)
            {
                throw new Exception("Bad x in routine betai");
            }

            if (x == 0.0 || x == 1.0)
                return x;

            if (a > SWITCH && b > SWITCH)
            {
                return betaiapprox(a, b, x);
            };

            bt = Math.Exp(gammln(a + b) - gammln(a) - gammln(b) + a * Math.Log(x) + b * Math.Log(1.0 - x));

            if (x < (a + 1.0) / (a + b + 2.0))
                return bt * betacf(a, b, x) / a;
            else
                return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
        }

        double betaiapprox(double a, double b, double x)
        {
            int j;

            double xu, t, sum, ans;

            double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);

            double lnmu = Math.Log(mu), lnmuc = Math.Log(1 - mu);

            t = Math.Sqrt(a * b / (Math.Pow((a + b), 2) * (a + b + 1.0)));

            if (x > a / (a + b))
            {
                if (x >= 1.0)
                    return 1.0;

                xu = Math.Min(1, Math.Max(mu + 10 * t, x + 5.0 * t));
            }
            else
            {
                if (x <= 0.0)
                    return 0.0;

                xu = Math.Max(0, Math.Min(mu - 10 * t, x - 5.0 * t));
            }

            sum = 0;

            for (j = 0; j < 18; j++)
            {
                t = x + (xu - x) * y[j];

                sum += w[j] * Math.Exp(a1 * (Math.Log(t) - lnmu) + b1 * (Math.Log(1 - t) - lnmuc));
            }

            ans = sum * (xu - x) * Math.Exp(a1 * lnmu - gammln(a) + b1 * lnmuc - gammln(b) + gammln(a + b));

            return ans > 0.0 ? 1.0 - ans : -ans;
        }

        /*
        #
        #  Inverce of incomplete Beta function
        #
        #
        */

        double invbetai(double p, double a, double b)
        {
            const double EPS = 1E-8; // 1.e - 8;

            double pp, t, u, err, x, al, h, w, afac, a1 = a - 1, b1 = b - 1;

            int j;

            if (p <= 0)
                return 0;
            else if (p >= 1)
                return 1;
            else if (a >= 1 && b >= 1)
            {
                pp = (p < 0.5) ? p : 1 - p;
                t = Math.Sqrt(-2 * Math.Log(pp));

                x = (2.30753 + t * 0.27061) / (1 + t * (0.99229 + t * 0.04481)) - t;

                if (p < 0.5)
                    x = -x;

                al = (Math.Pow(x, 2) - 3) / 6;

                h = 2 / (1 / (2 * a - 1) + 1 / (2 * b - 1));

                w = (x * Math.Sqrt(al + h) / h) - (1 / (2 * b - 1) - 1 / (2 * a - 1)) * (al + 5 / 6 - 2 / (3 * h));

                x = a / (a + b * Math.Exp(2 * w));
            }
            else
            {
                double lna = Math.Log(a / (a + b)), lnb = Math.Log(b / (a + b));

                t = Math.Exp(a * lna) / a;

                u = Math.Exp(b * lnb) / b;

                w = t + u;

                if (p < t / w)
                    x = Math.Pow(a * w * p, 1 / a);
                else
                    x = 1 - Math.Pow(b * w * (1 - p), 1 / b);
            }

            afac = -gammln(a) - gammln(b) + gammln(a + b);

            for (j = 0; j < 10; j++)
            {
                if (x == 0 || x == 1)
                    return x;

                err = betai(a, b, x) - p;

                t = Math.Exp(a1 * Math.Log(x) + b1 * Math.Log(1 - x) + afac);

                u = err / t;

                x -= (t = u / (1 - 0.5 * Math.Min(1, u * (a1 / x - b1 / (1 - x)))));

                if (x <= 0)
                    x = 0.5 * (x + t);

                if (x >= 1)
                    x = 0.5 * (x + t + 1);

                if (Math.Abs(t) < EPS * x && j > 0)
                    break;
            }
            return x;
        }

        /*
        #
        # nu - degrees of freedom
        #
        #
        */

        public double invcdf_T(double p, int nu)
        {
            if (p <= 0 || p >= 1)
                throw new Exception("bad p in Studentdist");

            double x;

            //printf("%10.5f %10.5f\n", 2.*MIN(p, 1. - p), 0.5*(double)nu);

            x = invbetai(2 * Math.Min(p, 1 - p), 0.5 * nu, 0.5);

            //x = sig*Math.Sqrt(nu*(1. - x) / x);

            x = Math.Sqrt(nu * (1 - x) / x);

            //return (p >= 0.5 ? mu + x : mu - x);

            return (p >= 0.5 ? x : -x);
        }

    }

}
