using UnityEngine;

public class BesselSystem
{
    float[][] _bessel;
    float[,] _zeroes;
    float[,] _norms;

    public BesselSystem(int maxHarmonic, int rootsNum)
    {
        _bessel = new float[maxHarmonic][];
        _zeroes = new float[maxHarmonic, rootsNum];
        _norms = new float[maxHarmonic, rootsNum];

        Debug.Log($"Begin init: {Time.realtimeSinceStartup:F3}");
        int N_0 = 20;
        int rootIter = 32; // prec = pi / 2 / 2^rootIter / 2 = pi / 2^(rootIter + 2)
        float prevRoot = 2.4048f + 0.1f - Mathf.PI / 2;
        for (int harmonic = 0; harmonic < maxHarmonic; harmonic++)
        {
            // space between first roots is < Mathf.PI / 2, J_n > 0 from 0 to first root => J_n < 0 on init
            float firstRoot = prevRoot + Mathf.PI / 2;
            for (int i = 0; i < rootIter; i++)
            {
                float mid = (prevRoot + firstRoot) * 0.5f;
                if (CalcBesselIntegrative(harmonic, mid) > 0)
                    prevRoot = mid;
                else
                    firstRoot = mid;
            }
            firstRoot = (prevRoot + firstRoot) * 0.5f;
            _zeroes[harmonic, 0] = firstRoot;
            prevRoot = firstRoot;

            float leftSign = -1f;
            for (int root = 1; root < rootsNum; root++)
            {
                // space between roots is > Mathf.PI for n > 0
                float nextRoot = firstRoot + Mathf.PI; // -0.0, 0.043, 0.141, 0.240, 0.335, 0.426
                while (leftSign * CalcBesselIntegrative(harmonic, nextRoot) > 0)
                    nextRoot += Mathf.PI * 0.5f;
                // find root
                for (int i = 0; i < rootIter; i++)
                {
                    float mid = (firstRoot + nextRoot) * 0.5f;
                    if (leftSign * CalcBesselIntegrative(harmonic, mid) > 0)
                        firstRoot = mid;
                    else
                        nextRoot = mid;
                }
                firstRoot = (firstRoot + nextRoot) * 0.5f;
                // write root
                _zeroes[harmonic, root] = firstRoot;
                firstRoot = nextRoot;
                leftSign *= -1f;
            }
        }
        Debug.Log($"Root inited: {Time.realtimeSinceStartup:F3}");

        for (int harmonic = 0; harmonic < maxHarmonic; harmonic++)
        {
            int size = N_0 * (harmonic + rootsNum) + 1; // include x = 0 and x = maxRoot
            _bessel[harmonic] = new float[size];
            float step = _zeroes[harmonic, rootsNum - 1] / size;
            float r = 0f; // step * 0.5f; => bad interpolation
                          // calc bessel to last root
            for (int i = 0; i < size; i++)
            {
                _bessel[harmonic][i] = CalcBesselIntegrative(harmonic, r);
                r += step;
            }
            // non-recursively calc norm for roots, ||J_n(mu_{n,k}/x_0 * x)||_[0, x_0] = x_0^2/2 (J_{n+1}(mu_{n,k}))^2, x_0 = 1
            // - http://www.lib.unn.ru/students/src/bessel.pdf
            for (int root = 0; root < rootsNum; root++)
            {
                float j_n1 = CalcBesselIntegrative(harmonic + 1, _zeroes[harmonic, root]);
                _norms[harmonic, root] = (2 * Mathf.PI) * 0.5f * j_n1 * j_n1;
                if (harmonic > 0)
                    _norms[harmonic, root] *= 0.5f;
            }
        }

        Debug.Log($"End init: {Time.realtimeSinceStartup:F3}");
        //return CalcBessel(index, r);
        //return CalcBesselIntegrative(index, _zeroes[index, zeroNum] * r01);
    }

    public float GetBesselNorm(int harmonic, int rootIndex)
    {
        return _norms[harmonic, rootIndex];
    }

    /** r from 0 to 1 */
    public float Bessel(int index, int zeroNum, float r01)
    {
        return Bessel(index, _zeroes[index, zeroNum] * r01);
    }
    public float Bessel(int index, float r)
    {
        if (index >= _zeroes.GetLength(0))
            return 0;
        float maxR = _zeroes[index, _zeroes.GetLength(1) - 1];
        if (r > maxR)
            return 0;

        float step = maxR / (_bessel[index].Length - 1);
        int highIndex = Mathf.CeilToInt(r / step);
        int lowIndex = highIndex - 1;
        if (lowIndex < 0)
            lowIndex = 0;
        float topDist = highIndex - r / step;
        float botDist = 1 - topDist;
        return _bessel[index][lowIndex] * topDist + _bessel[index][highIndex] * botDist;
    }

    static float CalcBessel(float index, float r)
    {
        float res;

        if (r > 200)
        {
            // use approx for large r and index
            float phase = (index + 0.5f) * Mathf.PI * 0.5f;
            res = Mathf.Sqrt(2f / (Mathf.PI * r)) * Mathf.Cos(r - phase);
        }
        else if (r > 20)
        {
            // sqrt(2 / pi r) [ cos (r - pi index / 2 - pi / 4) sum (from m = 0) (-1)^m (index, 2m) / (2m)^(2m)
            //     - sin (r - pi index / 2 - pi / 4) sum (from m = 0) (-1)^m (index, 2m + 1) / (2r)^(2m + 1) ]
            // (index, m) = Г (index + m + 0.5) / Г (m + 1) / Г (index - m + 0.5) = 1 / m! / prod (from k = -m + 1 to m) (index + k + 0.5)
            float cosCoef = 0;
            for (int m = 0; m <= 6; m++)
            {
                float c = 1;
                int m2 = 2 * m;
                for (int k = 2; k <= m2; k++)
                    c /= k;
                for (int k = -m2 + 1; k <= m2; k++)
                    c /= index + k + 0.5f;
                cosCoef += (1 - 2 * (m % 2)) * c / Mathf.Pow(m2, m2);
            }
            float sinCoef = 0;
            for (int m = 0; m <= 6; m++)
            {
                float c = 1;
                int m2 = 2 * m + 1;
                for (int k = 2; k <= m2; k++)
                    c /= k;
                for (int k = -m2 + 1; k <= m2; k++)
                    c /= index + k + 0.5f;
                sinCoef += (1 - 2 * (m % 2)) * c / Mathf.Pow(2 * r, m2);
            }
            float phase = r - (index + 0.5f) * Mathf.PI * 0.5f;
            res = cosCoef * Mathf.Cos(phase) - sinCoef * Mathf.Sin(phase);
            res = Mathf.Sqrt(2f / (Mathf.PI * r)) * res;
        }
        else
        {
            // sum from m = 0 to inf (-1)^m / (m! Г(m + index + 1)) (r/2)^(2m + index)

            float member = 1 / CalcGamma(index + 1) * Mathf.Pow(r * 0.5f, index);
            float rd = (r * 0.5f) * (r * 0.5f);
            res = member;
            for (int m = 1; m <= 30; m++)
            {
                member *= -1f / m / (m + index) * rd;
                res += member;
            }
        }

        return res;
    }

    static float CalcBesselIntegrative(int index, float r)
    {
        // 1 / 2pi int (from 0 to 2pi) e^(ir sin theta - in theta) d theta
        float res = 0;

        int N_0 = 20; // points on one period
        int N = N_0 * (index + Mathf.FloorToInt(r / Mathf.PI) + 1);
        float step = 2 * Mathf.PI / N;
        float c = 1f / N;
        for (int i = 0; i < N; i++)
        {
            float theta = step * i;
            res += Mathf.Cos(r * Mathf.Sin(theta) - index * theta) * c;
        }

        return res;
    }

    // 1 => 1
    // 2 => 1
    // 3 => 2
    // 4 => 6
    static float CalcGamma(float t)
    {
        // calc gamma (1, 2), then * m, because m Г(m) = Г(m+1)
        float m = t - Mathf.Floor(t) + 1;

        float res = 0;
        float x = 0;
        for (int k = 1; k <= 50; k++)
        {
            float nextX = 0.01f * k * k;
            float middleX = (x + nextX) * 0.5f;
            float deltaX = nextX - x;
            res += Mathf.Pow(middleX, m - 1) * Mathf.Exp(-middleX) * deltaX;
            x = nextX;
        }

        while (m < t - 0.1f)
        {
            res *= m;
            m += 1f;
        }

        //Debug.Log(t + " " + m + " " + res);
        return res;
    }
}