#ifndef BECPP_GRID_H
#define BECPP_GRID_H

#include <utility>
#include <vector>

using doubleArray_t = std::vector<std::vector<double>>;

class Grid
{
protected:
    virtual void constructGridParams() = 0;
    virtual void constructGrids() = 0;
    virtual void fftshift() = 0;
    virtual ~Grid() = default;
};

class Grid1D : public Grid
{
private:
    void constructGridParams() override;
    void constructGrids() override;
    void fftshift() override;

    const unsigned int m_xPoints{};
    double m_xGridSpacing{};
    double m_xFourierGridSpacing{};
    double m_xLength{};

    std::vector<double> m_xMesh{};
    std::vector<double> m_xFourierMesh{};
    std::vector<double> m_wavenumber{};

public:
    Grid1D(unsigned int nx, double dx);
    ~Grid1D() override = default;

    [[nodiscard]] unsigned int xPoints() const;
    [[nodiscard]] double xGridSpacing() const;
    [[nodiscard]] double xFourierGridSpacing() const;
    [[nodiscard]] double xLength() const;
    [[nodiscard]] std::vector<double> wavenumber() const;

    friend class Wavefunction;
};

class Grid2D : public Grid
{
private:
    void constructGridParams() override;
    void constructGrids() override;
    void fftshift() override;

    const unsigned int m_xPoints{};
    const unsigned int m_yPoints{};
    double m_xGridSpacing{};
    double m_yGridSpacing{};
    double m_gridSpacingProduct{};
    double m_xFourierGridSpacing{};
    double m_yFourierGridSpacing{};
    double m_xLength{};
    double m_yLength{};

    doubleArray_t m_xMesh{};
    doubleArray_t m_yMesh{};
    doubleArray_t m_xFourierMesh{};
    doubleArray_t m_yFourierMesh{};
    doubleArray_t m_wavenumber{};

public:
    Grid2D(unsigned int nx, unsigned int ny, double dx, double dy);
    ~Grid2D() override = default;

    [[nodiscard]] unsigned int xPoints() const;
    [[nodiscard]] unsigned int yPoints() const;
    [[nodiscard]] double xGridSpacing() const;
    [[nodiscard]] double yGridSpacing() const;
    [[nodiscard]] double gridSpacingProduct() const;
    [[nodiscard]] double xFourierGridSpacing() const;
    [[nodiscard]] double yFourierGridSpacing() const;
    [[nodiscard]] double xLength() const;
    [[nodiscard]] double yLength() const;
    [[nodiscard]] doubleArray_t wavenumber() const;

    friend class Wavefunction;
};

#endif //BECPP_GRID_H
