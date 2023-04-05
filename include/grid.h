#ifndef BECPP_GRID_H
#define BECPP_GRID_H

#include <tuple>
#include <utility>
#include <vector>

using vector2D_t = std::vector<std::vector<double>>;

class Grid
{
protected:
    virtual void constructGridParams() = 0;
    virtual void constructMesh() = 0;
    virtual void fftshift() = 0;
    virtual ~Grid() = default;
};

struct Mesh1D
{
    std::vector<double> m_xMesh{};
    std::vector<double> m_xFourierMesh{};
    std::vector<double> m_wavenumber{};
};

class Grid1D : public Grid
{
private:
    void constructGridParams() override;
    void constructMesh() override;
    void fftshift() override;

    const unsigned int m_gridPoints{};
    double m_gridSpacing{};
    double m_fourierGridSpacing{};
    double m_length{};
    Mesh1D m_mesh{};

public:
    Grid1D(unsigned int nx, double dx);
    ~Grid1D() override = default;

    [[nodiscard]] unsigned int shape() const;
    [[nodiscard]] double gridSpacing() const;
    [[nodiscard]] double fourierGridSpacing() const;
    [[nodiscard]] double gridLength() const;
    [[nodiscard]] std::vector<double> wavenumber() const;

    friend class Wavefunction;
};

struct Mesh2D
{
    vector2D_t xMesh{};
    vector2D_t yMesh{};
    vector2D_t xFourierMesh{};
    vector2D_t yFourierMesh{};
    vector2D_t wavenumber{};
};

class Grid2D : public Grid
{
private:
    void constructGridParams() override;
    void constructMesh() override;
    void fftshift() override;

    const std::tuple<unsigned int, unsigned int> m_gridPoints{};
    const std::tuple<double, double> m_gridSpacing{};
    std::tuple<double, double> m_fourierGridSpacing{};
    std::tuple<double, double> m_gridLength{};
    Mesh2D m_mesh{};

public:
    Grid2D(std::tuple<unsigned int, unsigned int> points,
           std::tuple<double, double> gridSpacing);
    ~Grid2D() override = default;

    [[nodiscard]] std::tuple<unsigned int, unsigned int> ndim() const;
    [[nodiscard]] std::tuple<double, double> gridSpacing() const;
    [[nodiscard]] std::tuple<double, double> fourierGridSpacing() const;
    [[nodiscard]] std::tuple<double, double> gridLength() const;
    [[nodiscard]] vector2D_t wavenumber() const;

    friend class Wavefunction;
};

class Grid3D : public Grid
{
private:
    void constructGridParams() override;
    void constructMesh() override;
    void fftshift() override;

    const unsigned int m_xPoints{};
    const unsigned int m_yPoints{};
    const unsigned int m_zPoints{};
    double m_xGridSpacing{};
    double m_yGridSpacing{};
    double m_zGridSpacing{};
    double m_gridSpacingProduct{};
    double m_xFourierGridSpacing{};
    double m_yFourierGridSpacing{};
    double m_zFourierGridSpacing{};
    double m_xLength{};
    double m_yLength{};
    double m_zLength{};

    std::vector<vector2D_t> m_xMesh{};
    std::vector<vector2D_t> m_yMesh{};
    std::vector<vector2D_t> m_zMesh{};
    std::vector<vector2D_t> m_xFourierMesh{};
    std::vector<vector2D_t> m_yFourierMesh{};
    std::vector<vector2D_t> m_zFourierMesh{};
    std::vector<vector2D_t> m_wavenumber{};

public:
    Grid3D(std::tuple<unsigned int, unsigned int, unsigned int> points,
           std::tuple<double, double, double> gridSpacing);
    ~Grid3D() override = default;

    [[nodiscard]] unsigned int xPoints() const;
    [[nodiscard]] unsigned int yPoints() const;
    [[nodiscard]] unsigned int zPoints() const;
    [[nodiscard]] double xGridSpacing() const;
    [[nodiscard]] double yGridSpacing() const;
    [[nodiscard]] double zGridSpacing() const;
    [[nodiscard]] double gridSpacingProduct() const;
    [[nodiscard]] double xFourierGridSpacing() const;
    [[nodiscard]] double yFourierGridSpacing() const;
    [[nodiscard]] double zFourierGridSpacing() const;
    [[nodiscard]] double xLength() const;
    [[nodiscard]] double yLength() const;
    [[nodiscard]] double zLength() const;
    [[nodiscard]] std::vector<vector2D_t> wavenumber() const;

    friend class Wavefunction;
};

#endif//BECPP_GRID_H
