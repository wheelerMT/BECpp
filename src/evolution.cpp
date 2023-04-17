#include "evolution.h"

void fourierStep(Wavefunction1D& wfn, const Parameters& params)
{
    for (int i = 0; i < wfn.grid().shape(); ++i)
    {
        wfn.component()[i] *=
                exp(-0.25 * I * params.timeStep * wfn.grid().wavenumber()[i]);
    }
}

void fourierStep(Wavefunction2D& wfn, const Parameters& params)
{
    auto [xPoints, yPoints] = wfn.grid().shape();
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            auto index = j + i * yPoints;
            wfn.component()[index] *= exp(-0.25 * I * params.timeStep *
                                          wfn.grid().wavenumber()[index]);
        }
    }
}

void fourierStep(Wavefunction3D& wfn, const Parameters& params)
{
    auto [xPoints, yPoints, zPoints] = wfn.grid().shape();
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int k = 0; k < zPoints; ++k)
            {
                auto index = k + zPoints * (j + i * yPoints);
                wfn.component()[index] *= exp(-0.25 * I * params.timeStep *
                                              wfn.grid().wavenumber()[index]);
            }
        }
    }
}

void interactionStep(Wavefunction1D& wfn, const Parameters& params)
{
    for (int i = 0; i < wfn.grid().shape(); ++i)
    {
        wfn.component()[i] *=
                exp(-I * params.timeStep *
                    (params.trap[i] +
                     params.intStrength * pow(abs(wfn.component()[i]), 2)));
    }
}

void interactionStep(Wavefunction2D& wfn, const Parameters& params)
{
    auto [xPoints, yPoints] = wfn.grid().shape();
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            auto index = j + i * yPoints;
            wfn.component()[index] *= exp(
                    -I * params.timeStep *
                    (params.trap[index] +
                     params.intStrength * pow(abs(wfn.component()[index]), 2)));
        }
    }
}

void interactionStep(Wavefunction3D& wfn, const Parameters& params)
{
    auto [xPoints, yPoints, zPoints] = wfn.grid().shape();
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int k = 0; k < zPoints; ++k)
            {
                auto index = k + zPoints * (j + i * yPoints);
                wfn.component()[index] *=
                        exp(-I * params.timeStep *
                            (params.trap[index] +
                             params.intStrength *
                                     pow(abs(wfn.component()[index]), 2)));
            }
        }
    }
}

double calculateAtomNum(const Wavefunction1D& wfn)
{
    double atomNumber{};
    std::vector<double> dens = wfn.density();
    for (int i = 0; i < wfn.grid().shape(); ++i)
    {
        atomNumber += std::pow(std::abs(dens[i]), 2) * wfn.grid().gridSpacing();
    }

    return atomNumber;
}

double calculateAtomNum(const Wavefunction2D& wfn)
{
    double atomNumber{};
    auto [xPoints, yPoints] = wfn.grid().shape();
    auto [xGridSpacing, yGridSpacing] = wfn.grid().gridSpacing();
    std::vector<double> dens = wfn.density();

    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            auto index = j + i * yPoints;
            atomNumber += std::pow(std::abs(dens[index]), 2) * xGridSpacing *
                          yGridSpacing;
        }
    }

    return atomNumber;
}

double calculateAtomNum(const Wavefunction3D& wfn)
{
    double atomNumber{};
    auto [xPoints, yPoints, zPoints] = wfn.grid().shape();
    auto [xGridSpacing, yGridSpacing, zGridSpacing] = wfn.grid().gridSpacing();
    std::vector<double> dens = wfn.density();

    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int k = 0; k < zPoints; ++k)
            {
                auto index = k + zPoints * (j + yPoints * i);
                atomNumber += std::pow(std::abs(dens[index]), 2) *
                              xGridSpacing * yGridSpacing * zGridSpacing;
            }
        }
    }

    return atomNumber;
}

void renormaliseAtomNum(Wavefunction1D& wfn)
{
    double currentAtomNum = calculateAtomNum(wfn);

    for (int i = 0; i < wfn.grid().shape(); ++i)
    {
        wfn.component()[i] *= sqrt(wfn.atomNumber()) / sqrt(currentAtomNum);
    }
}

void renormaliseAtomNum(Wavefunction2D& wfn)
{
    double currentAtomNum = calculateAtomNum(wfn);
    auto [xPoints, yPoints] = wfn.grid().shape();

    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            auto index = j + i * yPoints;
            wfn.component()[index] *=
                    sqrt(wfn.atomNumber()) / sqrt(currentAtomNum);
        }
    }
}

void renormaliseAtomNum(Wavefunction3D& wfn)
{
    double currentAtomNum = calculateAtomNum(wfn);
    auto [xPoints, yPoints, zPoints] = wfn.grid().shape();

    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int k = 0; k < zPoints; ++k)
            {
                auto index = k + zPoints * (j + i * yPoints);
                wfn.component()[index] *=
                        sqrt(wfn.atomNumber()) / sqrt(currentAtomNum);
            }
        }
    }
}
