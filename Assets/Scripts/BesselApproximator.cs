using System;
using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEngine;
using UnityEngine.UI;

public class BesselApproximator : MonoBehaviour
{
    public Texture2D Original;
    public Vector2Int OutputSize = new Vector2Int(400, 400);
    Texture2D _display;
    float _prevTime = 0;

    BesselSystem _besselSystem; // TODO: static system
    IRadialImageData _imageData;

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
    }

    void Update()
    {
        if (Time.time - _prevTime < 1f)
            return;

        _prevTime = Time.time;
        Draw(_display, _imageData);
    }

    void Draw(Texture2D texture, IRadialImageData data)
    {
        Debug.Log($"Begin drawing: {Time.realtimeSinceStartup:F3}");
        Color[] cols = texture.GetPixels();

        float stepX = 2 / (float)texture.width;
        float stepY = 2 / (float)texture.height;
        float y = -1f + stepY * 0.5f;
        for (int iy = 0; iy < texture.height; ++iy)
        {
            float x = -1f + stepX * 0.5f;
            for (int ix = 0; ix < texture.width; ++ix)
            {
                float r = Mathf.Sqrt(x * x + y * y);
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

    IRadialImageData Convert(Texture2D original, BesselSystem besselSystem, RadialImagePrecision prec, BesselConvertingMode mode)
    {
        if (mode == BesselConvertingMode.RGB)
        {
            Debug.Log($"Begin converting: {Time.realtimeSinceStartup:F3}");
            BesselSeriesData redSeriesData = new BesselSeriesData(prec.HarmonicCount, prec.RootsCount);
            BesselSeriesData greenSeriesData = new BesselSeriesData(prec.HarmonicCount, prec.RootsCount);
            BesselSeriesData blueSeriesData = new BesselSeriesData(prec.HarmonicCount, prec.RootsCount);
            ConvertingData jobR = new ConvertingData(ConvertingData.ConvertingChannel.Red, original, redSeriesData, prec, besselSystem);
            ConvertingData jobG = new ConvertingData(ConvertingData.ConvertingChannel.Green, original, greenSeriesData, prec, besselSystem);
            ConvertingData jobB = new ConvertingData(ConvertingData.ConvertingChannel.Blue, original, blueSeriesData, prec, besselSystem);
            Parallel.For(0, prec.HarmonicCount * prec.RootsCount, (i) => jobR.Execute(i));
            Parallel.For(0, prec.HarmonicCount * prec.RootsCount, (i) => jobG.Execute(i));
            Parallel.For(0, prec.HarmonicCount * prec.RootsCount, (i) => jobB.Execute(i));
            Debug.Log($"End converting: {Time.realtimeSinceStartup:F3}");

            return new RGBRadialImageData(redSeriesData, greenSeriesData, blueSeriesData, besselSystem);
        }
        else
        {
            Debug.Log($"Begin converting: {Time.realtimeSinceStartup:F3}");
            BesselSeriesData seriesData = new BesselSeriesData(prec.HarmonicCount, prec.RootsCount);
            ConvertingData job = new ConvertingData(ConvertingData.ConvertingChannel.Red, original, seriesData, prec, besselSystem);
            Parallel.For(0, prec.HarmonicCount * prec.RootsCount, (i) => job.Execute(i));
            Debug.Log($"End converting: {Time.realtimeSinceStartup:F3}");

            return new GrayscaleRadialImageData(seriesData, _besselSystem);
        }
    }

    class ConvertingData
    {
        BesselSystem _besselSystem;
        RadialImagePrecision _prec;
        BesselSeriesData _seriesData;
        Color[] _pixels;
        int width, height;
        ConvertingChannel _channel;

        public ConvertingData(ConvertingChannel channel, Texture2D original, BesselSeriesData seriesData, RadialImagePrecision prec, BesselSystem besselSystem)
        {
            _seriesData = seriesData;
            _prec = prec;
            _besselSystem = besselSystem;
            width = original.width;
            height = original.height;
            _pixels = original.GetPixels();
            _channel = channel;
        }

        public void Execute(int i)
        {
            int n = i / _prec.RootsCount;
            int k = i % _prec.RootsCount;
            float cosCoef = 0f;
            float sinCoef = 0f;
            // 2D integral r J cos
            int rPoints = _prec.IntegralPointsR(n, k, width, height); // 400x400 => 400 * 1.4
            float dr = 1 / (float)rPoints;
            for (int rIndex = 0; rIndex < rPoints; rIndex++)
            {
                float r01 = dr * (rIndex + 0.5f);
                float bessel = _besselSystem.Bessel(n, k, r01);
                bessel /= _besselSystem.GetBesselNorm(n, k);

                int phiPoints = _prec.IntegralPointsPhi(n, k, width, height, r01); // 400x400 => 1 + 8 + 14 + ... + (1 + 6.28 * 400 * 1.4)
                // 400x400 => ~400 + 2 * pi * sqrt(2) * (400 * (400 - 1)) / 2 ~ pi * sqrt(2) * 400 * 400 > 400 * 400 => can use blurred images
                float dphi = 2f * Mathf.PI / (float)phiPoints;
                float impact = r01 * dr * dphi; // [r] dr dphi = weight
                for (int phiIndex = 0; phiIndex < phiPoints; phiIndex++)
                {
                    float phi = dphi * phiIndex;
                    int pixelX = Mathf.FloorToInt(width * (0.5f + 0.5f * r01 * Mathf.Cos(phi)));
                    int pixelY = Mathf.FloorToInt(height * (0.5f + 0.5f * r01 * Mathf.Sin(phi)));
                    int pixelIndex = pixelY * height + pixelX;
                    float color;
                    switch (_channel)
                    {
                        case ConvertingChannel.Red: color = _pixels[pixelIndex].r; break;
                        case ConvertingChannel.Green: color = _pixels[pixelIndex].g; break;
                        case ConvertingChannel.Blue: color = _pixels[pixelIndex].b; break;
                        default: color = 0f; break;
                    }
                    color = 2f * color - 1f;
                    cosCoef += color * bessel * Mathf.Cos(n * phi) * impact;
                    sinCoef += color * bessel * Mathf.Sin(n * phi) * impact;
                }
            }
            _seriesData.CosFourierCoefficients[n, k] = cosCoef;
            _seriesData.SinFourierCoefficients[n, k] = sinCoef;
        }

        public enum ConvertingChannel
        {
            Red, Green, Blue
        }
    }
}
