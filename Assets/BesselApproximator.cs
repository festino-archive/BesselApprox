using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class BesselApproximator : MonoBehaviour
{
    public Texture2D Original;
    public Vector2Int OutputSize = new Vector2Int(400, 400);
    Texture2D _display;
    float _prevTime = 0;

    BesselSystem _besselSystem; // TODO: static system
    RadialImageData _imageData;

    void Start()
    {
        //https://docs.unity3d.com/ScriptReference/Texture2D.SetPixels.html
        RawImage rend = GetComponent<RawImage>();

        // duplicate the original texture and assign to the material
        _display = new Texture2D(OutputSize.x, OutputSize.y);
        rend.texture = _display;

        _besselSystem = new BesselSystem(30, 30); // TODO: static system
        RadialImagePrecision prec = new RadialImagePrecision(20, 20);
        _imageData = Convert(Original, _besselSystem, prec, BesselConvertingMode.GRAYSCALE);
        /*BesselSeriesData seriesData = new BesselSeriesData(prec.HarmonicCount, prec.RootsCount);
        seriesData.CosFourierCoefficients[0, 4] = 0.5f;
        seriesData.CosFourierCoefficients[1, 7] = -0.4f;
        seriesData.CosFourierCoefficients[3, 9] = 1f;
        seriesData.SinFourierCoefficients[7, 9] = 1f;
        imageData = new GrayscaleRadialImageData(seriesData);*/
    }

    void Update()
    {
        if (Time.time - _prevTime < 1f)
            return;

        _prevTime = Time.time;
        Draw(_display, _imageData);
    }

    RadialImageData Convert(Texture2D original, BesselSystem besselSystem, RadialImagePrecision prec, BesselConvertingMode mode)
    {
        if (mode == BesselConvertingMode.RBG)
        {
            return new RGBRadialImageData();
        }
        else
        {
            Debug.Log($"Begin converting: {Time.realtimeSinceStartup:F3}");
            BesselSeriesData seriesData = new BesselSeriesData(prec.HarmonicCount, prec.RootsCount);
            for (int n = 0; n < prec.HarmonicCount; n++)
            {
                for (int k = 0; k < prec.RootsCount; k++)
                {
                    float cosCoef = 0;
                    float sinCoef = 0;
                    // 2D integral r J cos
                    int rPoints = prec.IntegralPointsR(n, k, original.width, original.height);
                    float dr = 1 / (float)rPoints;
                    for (int rIndex = 0; rIndex < rPoints; rIndex++)
                    {
                        float r01 = dr * (rIndex + 0.5f);
                        float bessel = besselSystem.Bessel(n, k, r01);
                        bessel /= _besselSystem.GetBesselNorm(n, k);

                        int phiPoints = prec.IntegralPointsPhi(n, k, original.width, original.height, r01);
                        float dphi = 2 * Mathf.PI / (float)phiPoints;
                        float impact = r01 * dr * dphi; // [r] dr dphi = weight
                        for (int phiIndex = 0; phiIndex < phiPoints; phiIndex++)
                        {
                            float phi = dphi * phiIndex;
                            int pixelX = Mathf.FloorToInt(original.width * (0.5f + 0.5f * r01 * Mathf.Cos(phi)));
                            int pixelY = Mathf.FloorToInt(original.height * (0.5f + 0.5f * r01 * Mathf.Sin(phi)));
                            float color = 2 * original.GetPixel(pixelX, pixelY).r - 1f;
                            cosCoef += color * bessel * Mathf.Cos(n * phi) * impact;
                            sinCoef += color * bessel * Mathf.Sin(n * phi) * impact;
                        }
                    }
                    seriesData.CosFourierCoefficients[n, k] = cosCoef;
                    seriesData.SinFourierCoefficients[n, k] = sinCoef;
                }
            }
            Debug.Log($"End converting: {Time.realtimeSinceStartup:F3}");
            return new GrayscaleRadialImageData(seriesData, _besselSystem);
        }
    }

    void Draw(Texture2D texture, RadialImageData data)
    {
        Debug.Log($"Begin drawing: {Time.realtimeSinceStartup:F3}");
        Color[] cols = texture.GetPixels();
        //for (int i = 0; i < cols.Length; ++i)
        //    cols[i] = new Color((i % 256) / 255.0f, (i / 256) % 256 / 255.0f, (i / 256 / 256) % 256 / 255.0f);

        float stepX = 2 / (float)texture.width;
        float stepY = 2 / (float)texture.height;
        float y = stepY * 0.5f - stepY * texture.height * 0.5f;
        for (int iy = 0; iy < texture.height; ++iy)
        {
            float x = stepX * 0.5f - stepX * texture.width * 0.5f;
            for (int ix = 0; ix < texture.width; ++ix)
            {
                //cols[iy * texture.height + ix] = new Color(ix / (float)texture.width, iy / (float)texture.height, 0);
                //cols[iy * texture.height + ix] = data.GetColor(2 * ix / (float)texture.width - 1, 2 * iy / (float)texture.height - 1);
                float r = Mathf.Sqrt(x * x + y * y);// / Mathf.Sqrt(2);
                float phi = Mathf.Atan2(y, x);
                cols[iy * texture.height + ix] = data.GetColor(r, phi);
                x += stepX;
            }
            y += stepY;
        }
        texture.SetPixels(cols);
        texture.Apply(true);
        Debug.Log($"End drawing: {Time.realtimeSinceStartup:F3}");
    }

    interface RadialImageData
    {
        public Color GetColor(float r, float phi);
    }
    class RGBRadialImageData : RadialImageData
    {
        public Color GetColor(float r, float phi)
        {
            return new Color(r, 1-r, r);
        }
    }
    class GrayscaleRadialImageData : RadialImageData
    {
        BesselSystem _besselSystem;
        BesselSeriesData _seriesData; // one channel

        public GrayscaleRadialImageData(BesselSeriesData seriesData, BesselSystem system)
        {
            _seriesData = seriesData;
            _besselSystem = system;
        }

        public Color GetColor(float r01, float phi)
        {
            float res = 0;
            for (int harmonic = 0; harmonic < _seriesData.HarmonicCount; harmonic++)
            {
                float cos = Mathf.Cos(harmonic * phi);
                float sin = Mathf.Sin(harmonic * phi);
                for (int rootIndex = 0; rootIndex < _seriesData.RootsCount; rootIndex++)
                {
                    float bessel = _besselSystem.Bessel(harmonic, rootIndex, r01);
                    float member = _seriesData.CosFourierCoefficients[harmonic, rootIndex] * cos * bessel;
                    member += _seriesData.SinFourierCoefficients[harmonic, rootIndex] * sin * bessel;
                    res += member;
                }
            }
            res = 0.5f + 0.5f * res;
            return new Color(res, res, res);
        }
    }

    //Prec is: (Bessel calc settings), number of harmonics, number of Bessel_n roots, integral resolution
    struct RadialImagePrecision
    {
        public int HarmonicCount;
        public int RootsCount;

        public RadialImagePrecision(int harmonicCount, int rootsCount)
        {
            HarmonicCount = harmonicCount;
            RootsCount = rootsCount;
        }

        public int IntegralPointsR(int n, int k, int width, int height)
        {
            return Mathf.CeilToInt(Mathf.Sqrt(width * width + height * height));
        }

        public int IntegralPointsPhi(int n, int k, int width, int height, float r01)
        {
            return Mathf.CeilToInt(1 + 2 * Mathf.PI * r01 * Mathf.Sqrt(width * width + height * height));
        }
    }

    enum BesselConvertingMode
    {
        RBG, GRAYSCALE
    }

    struct BesselSeriesData
    {
        public int HarmonicCount { get => CosFourierCoefficients.GetLength(0); }
        public int RootsCount { get => CosFourierCoefficients.GetLength(1); }

        public readonly float[,] CosFourierCoefficients;
        public readonly float[,] SinFourierCoefficients;

        public BesselSeriesData(int harmCount, int rootCount)
        {
            CosFourierCoefficients = new float[harmCount, rootCount];
            SinFourierCoefficients = new float[harmCount, rootCount];
        }
    }

    class BesselSystem
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
                    float nextRoot = firstRoot + Mathf.PI; // -..., 0.043, 0.141, 0.240, 0.335, 0.426
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
                // (index, m) = Ã (index + m + 0.5) / Ã (m + 1) / Ã (index - m + 0.5) = 1 / m! / prod (from k = -m + 1 to m) (index + k + 0.5)
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
                // sum from m = 0 to inf (-1)^m / (m! Ã(m + index + 1)) (r/2)^(2m + index)

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
            // calc gamma (1, 2), then * m, because m Ã(m) = Ã(m+1)
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
}
