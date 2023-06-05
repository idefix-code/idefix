// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_PLANETARYSYSTEM_PLANET_HPP_
#define DATABLOCK_PLANETARYSYSTEM_PLANET_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "planetStructs.hpp"

// forward class declaration
class DataBlock;
class PlanetarySystem;


class Planet {
 public:
    Force m_force;
    Planet(int, Input &, DataBlock *, PlanetarySystem *);
    Planet(const Planet&);
    Planet& operator=(const Planet &p);
    //void Init(int &, Input &, DataBlock *, PlanetarySystem *);
    void displayPlanet() const;
    void RegisterInDump();
    void ShowConfig();
    real getMp() const;
    real getXp() const;
    real getYp() const;
    real getZp() const;
    real getVxp() const;
    real getVyp() const;
    real getVzp() const;
    bool getIsActive() const;
    int getIndex() const;
    void setMp(real qp);
    void setXp(real xp);
    void setYp(real yp);
    void setZp(real zp);
    void setVxp(real vxp);
    void setVyp(real vyp);
    void setVzp(real vzp);
    void updateMp(const real);
    void activatePlanet(const real);
    // refresh the force
    Point computeAccel(DataBlock&, bool&);
    void computeForce(DataBlock&, bool&);

 protected:
    friend class PlanetarySystem;
    DataBlock *data;
    PointSpeed state;

    real &m_xp;
    real &m_yp;
    real &m_zp;
    real &m_vxp;
    real &m_vyp;
    real &m_vzp;
    bool m_isActive;
    real m_qp;
    real m_qpIni;
    real m_tOffset;
    int m_ip;
    PlanetarySystem *pSys;
};

#endif // DATABLOCK_PLANETARYSYSTEM_PLANET_HPP_
