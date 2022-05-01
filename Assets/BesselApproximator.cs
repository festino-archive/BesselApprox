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
                    int rPoints = prec.IntegralPointsR(n, k, original.width, original.height); // 400x400 => 400 * 1.4
                    float dr = 1 / (float)rPoints;
                    for (int rIndex = 0; rIndex < rPoints; rIndex++)
                    {
                        float r01 = dr * (rIndex + 0.5f);
                        float bessel = besselSystem.Bessel(n, k, r01);
                        bessel /= _besselSystem.GetBesselNorm(n, k);

                        int phiPoints = prec.IntegralPointsPhi(n, k, original.width, original.height, r01); // 400x400 => 1 + 8 + 14 + ... + (1 + 6.28 * 400 * 1.4)
                        // 400x400 => ~400 + 2 * pi * sqrt(2) * (400 * (400 - 1)) / 2 ~ pi * sqrt(2) * 400 * 400 > 400 * 400 => can use blurred images
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
}
