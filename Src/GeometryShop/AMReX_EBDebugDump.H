
/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#ifndef _ALLREGULARSERVICE_H_
#define _ALLREGULARSERVICE_H_

#include "AMReX_GeometryService.H"

namespace amrex
{

/// A geometry service that is always regular.
/**
 */
  class AllRegularService: public GeometryService
  {
  public:
    AllRegularService();

    virtual ~AllRegularService();

    ///
    /**
       Return true if every cell in region is regular at the
       refinement described by dx.
    */
    virtual bool isRegular(const Box& a_region,
                           const Box& a_domain,
                           const RealVect& a_origin,
                           const Real& a_dx) const;

    ///
    /**
       Return true if every cell in region is covered at the
       refinement described by dx.
    */
    virtual bool isCovered(const Box& a_region,
                           const Box& a_domain,
                           const RealVect& a_origin,
                           const Real& a_dx) const;

    ///
    /**
     */
    virtual void fillGraph(BaseFab<int>&        a_regIrregCovered,
                           std::vector<IrregNode>&   a_nodes,
                           const Box&           a_validRegion,
                           const Box&           a_ghostRegion,
                           const Box& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;
  };
}
#endif
