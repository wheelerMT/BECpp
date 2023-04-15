#ifndef BECPP_GRID_H
#define BECPP_GRID_H

#include <tuple>
#include <utility>
#include <vector>

class Grid
{
protected:
    virtual void constructGridParams() = 0;
    virtual void constructMesh() = 0;
    virtual ~Grid() = default;
};

struct Mesh1D
{
    std::vector<double> xMesh{};
    std::vector<double> xFourierMesh{};
    std::vector<double> wavenumber{};
};

class Grid1D : public Grid
{
private:
    void constructGridParams() override;
    void constructMesh() override;

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
    std::vector<double> xMesh{};
    std::vector<double> yMesh{};
    std::vector<double> xFourierMesh{};
    std::vector<double> yFourierMesh{};
    std::vector<double> wavenumber{};
};

class Grid2D : public Grid
{
private:
    void constructGridParams() override;
    void constructMesh() override;

    const std::tuple<unsigned int, unsigned int> m_gridPoints{};
    const std::tuple<double, double> m_gridSpacing{};
    std::tuple<double, double> m_fourierGridSpacing{};
    std::tuple<double, double> m_gridLength{};
    Mesh2D m_mesh{};

public:
    Grid2D(std::tuple<unsigned int, unsigned int> points,
           std::tuple<double, double> gridSpacing);
    ~Grid2D() override = default;

    [[nodiscard]] std::tuple<unsigned int, unsigned int> shape() const;
    [[nodiscard]] std::tuple<double, double> gridSpacing() const;
    [[nodiscard]] std::tuple<double, double> fourierGridSpacing() const;
    [[nodiscard]] std::tuple<double, double> gridLength() const;
    [[nodiscard]] std::vector<double> wavenumber() const;

    friend class Wavefunction;
};

struct Mesh3D
{
    std::vector<double> xMesh{};
    std::vector<double> yMesh{};
    std::vector<double> zMesh{};
    std::vector<double> xFourierMesh{};
    std::vector<double> yFourierMesh{};
    std::vector<double> zFourierMesh{};
    std::vector<double> wavenumber{};
};

class Grid3D : public Grid
{
private:
    void constructGridParams() override;
    void constructMesh() override;

    const std::tuple<unsigned int, unsigned int, unsigned int>
            m_gridPoints{};
    const std::tuple<double, double, double> m_gridSpacing{};
    std::tuple<double, double, double> m_fourierGridSpacing{};
    std::tuple<double, double, double> m_gridLength{};
    Mesh3D m_mesh{};

public:
    Grid3D(std::tuple<unsigned int, unsigned int, unsigned int> points,
           std::tuple<double, double, double> gridSpacing);
    ~Grid3D() override = default;

    [[nodiscard]] std::tuple<unsigned int, unsigned int, unsigned int>
    shape() const;
    [[nodiscard]] std::tuple<double, double, double> gridSpacing() const;
    [[nodiscard]] std::tuple<double, double, double> fourierGridSpacing() const;
    [[nodiscard]] std::tuple<double, double, double> gridLength() const;
    [[nodiscard]] std::vector<double> wavenumber() const;

    friend class Wavefunction;
};

#endif//BECPP_GRID_H
