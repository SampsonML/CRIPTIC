#ifndef FLASHGG_H
#define FLASHGG_H

#include <mpi.h>
#include <hdf5.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "HDFIO_mpi.h"
#include "QUICKFLASH_HDF5.h" // for file IO in FLASH files
#include "GRID3D.h" // for uniform grid operations

// define float or double mode of FLASH GG ('general grid'; supposed to work for both UG and AMR)
// (note that float is the default because FLASH plt files are written in single precision)
// (if FLASH GG should use double precision,
//  both FLASH_GG_REAL and FLASH_GG_H5_REAL need to be defined accordingly by the user)
#ifndef FLASH_GG_REAL
#define FLASH_GG_REAL float
#define FLASH_GG_H5_REAL H5T_NATIVE_FLOAT
#endif

/**
 * FlashGG class
 * handels Uniform Grid and AMR (Paramesh) FLASH v3/4 files
 *
 * @author Christoph Federrath (christoph.federrath@anu.edu.au)
 * @version 2020
 *
 */

class FlashGG
{
    private:
    enum {X, Y, Z};
    std::string Inputfilename;
    std::string bounding_box_datasetname;
    std::string node_type_datasetname;
    char grid_type; // 'U' is UG, 'A' is AMR
    int NumBlocks, NumBlocksRep, NumDims, NBXY;
    std::vector<int> NB, NumBlocksIn, N;
    std::vector< std::vector< std::vector<double> > > BoundingBox;
    std::vector<int> NodeType;
    std::vector< std::vector<double> > MinMaxDomain, LBlock;
    std::vector< std::vector<double> > D;
    std::vector<double> Dmin;
    std::vector<double> Dmax;
    std::vector<double> L;
    bool Debug;
    // for pseudo blocks
    int NumBlocks_PB, NBXY_PB;
    std::vector<int> NB_PB, NumBlocksIn_PB;
    std::vector< std::vector<double> > LBlock_PB;
    std::vector< std::vector< std::vector<double> > > BoundingBox_PB;

    /**
      * Default constructor.
      */
    public: FlashGG(const std::string flashfile):
      Inputfilename(flashfile),
      bounding_box_datasetname("bounding box"),
      node_type_datasetname("node type"),
      NumBlocks(0), NumDims(0), NBXY(0),
      Debug(false)
    {
      // following calling sequence matters!
      if (Debug) { std::cout << "FlashGG: calling ReadNumBlocks..."<< std::endl;}
      this->ReadNumBlocks();
      if (Debug) { std::cout << "FlashGG: calling ReadNumCellsInBlock..."<< std::endl;}
      this->ReadNumCellsInBlock();
      if (Debug) { std::cout << "FlashGG: calling ReadBoundingBoxAndMinMaxDomain..."<< std::endl;}
      this->ReadBoundingBoxAndMinMaxDomain();
      if (Debug) { std::cout << "FlashGG: calling ReadNodeType..."<< std::endl;}
      this->ReadNodeType();
      if (Debug) { std::cout << "FlashGG: calling GetGridType..."<< std::endl;}
      grid_type = this->GetGridType();
      if (Debug) { std::cout << "FlashGG: Grid type is "<<grid_type<< std::endl;}
      // initialize pseudo blocks to be equal to actual blocks if UG
      if (grid_type == 'U') {
          // (user can call this from outside to set a requested number of cells per pseudo block)
          std::vector<int> ncells_pb = NB;
          if (Debug) { std::cout << "FlashGG: calling SetupPseudoBlocks..."<< std::endl;}
          this->SetupPseudoBlocks(ncells_pb);
      }
      if (Debug) {
          std::cout<<"FlashGG.h: FlashGG object created for file "<<flashfile<<"."<<std::endl;
          this->PrintInfo();
      }
    };

    public: void SetupPseudoBlocks(const std::vector<int> ncells_pb)
    {
        // divide whole domain in pseudo blocks
        // argument ncells_pb is the requested number of cells per pseudo block (per dimension)
        NB_PB.resize(NumDims);
        NumBlocksIn_PB.resize(NumDims);
        NumBlocks_PB = 1;
        for (int dim = 0; dim < NumDims; dim++) {
            NB_PB[dim] = ncells_pb[dim];
            // number of cells in pseudo blocks should always be integer multiples of total N cells in domain
            assert (N[dim] % NB_PB[dim] == 0);
            // define number of PBs along [X,Y,Z]
            NumBlocksIn_PB[dim] = N[dim] / NB_PB[dim];
            NumBlocks_PB *= NumBlocksIn_PB[dim];
        }
        // set up bounding box, etc for pseudo blocks
        NBXY_PB = NB_PB[X]*NB_PB[Y];
        BoundingBox_PB.resize(NumBlocks_PB);
        LBlock_PB.resize(NumBlocks_PB);
        for (int block = 0; block < NumBlocks_PB; block++) {
            LBlock_PB[block].resize(NumDims);
            BoundingBox_PB[block].resize(NumDims);
            for (int dim = 0; dim < NumDims; dim++) {
                LBlock_PB[block][dim] = L[dim] / NumBlocksIn_PB[dim];
                BoundingBox_PB[block][dim].resize(2);
            }
            // assume 3D blocks here
            int kmodb = block % (NumBlocksIn_PB[X]*NumBlocksIn_PB[Y]);
            int kb = block / (NumBlocksIn_PB[X]*NumBlocksIn_PB[Y]);
            int jb = kmodb / NumBlocksIn_PB[X];
            int ib = kmodb % NumBlocksIn_PB[X];
            BoundingBox_PB[block][X][0] = MinMaxDomain[X][0] + ib*LBlock_PB[block][X];
            BoundingBox_PB[block][X][1] = BoundingBox_PB[block][X][0] + LBlock_PB[block][X];
            BoundingBox_PB[block][Y][0] = MinMaxDomain[Y][0] + jb*LBlock_PB[block][Y];
            BoundingBox_PB[block][Y][1] = BoundingBox_PB[block][Y][0] + LBlock_PB[block][Y];
            BoundingBox_PB[block][Z][0] = MinMaxDomain[Z][0] + kb*LBlock_PB[block][Z];
            BoundingBox_PB[block][Z][1] = BoundingBox_PB[block][Z][0] + LBlock_PB[block][Z];
        }
    };

    /// PrintInfo (overloaded)
    public: void PrintInfo(void)
    {
        this->PrintInfo(false);
    }
    /// PrintInfo
    public: void PrintInfo(bool by_block)
    {
        std::cout<<"FlashGG: Number of dimensions (NumDims) = "<<NumDims<<std::endl;
        std::cout<<"FlashGG: Total number of blocks (NumBlocks) = "<<NumBlocks<<std::endl;
        if (grid_type == 'A') {
            std::vector<int> LeafBlocks = this->GetLeafBlocks();
            std::cout<<"FlashGG: Number of leaf blocks = "<<LeafBlocks.size()<<std::endl;
            std::cout<<"FlashGG: Max effective grid resolution: "<<N[X]<<" "<<N[Y]<<" "<<N[Z]<<std::endl;
        }
        if (grid_type == 'U') {
            std::cout<<"FlashGG: Total grid resolution (N) = "<<N[X]<<" "<<N[Y]<<" "<<N[Z]<<std::endl;
        }
        std::cout<<"FlashGG: Number of cells in block (NB) = "<<NB[X]<<" "<<NB[Y]<<" "<<NB[Z]<<std::endl;
        if (grid_type == 'U') {
            std::cout<<"FlashGG: Number of cells in pseudo block (NB_PB) = "<<NB_PB[X]<<" "<<NB_PB[Y]<<" "<<NB_PB[Z]<<std::endl;
            std::cout<<"FlashGG: Number of pseudo blocks in x,y,z (NumBlocksIn_PB) = "<<NumBlocksIn_PB[X]<<" "<<NumBlocksIn_PB[Y]<<" "<<NumBlocksIn_PB[Z]<<std::endl;
            std::cout<<"FlashGG: Total number of pseudo blocks (NumBlocks_PB) = "<<NumBlocks_PB<<std::endl;
        }
        std::cout<<"FlashGG: Min domain = "<<MinMaxDomain[X][0]<<" "<<MinMaxDomain[Y][0]<<" "<<MinMaxDomain[Z][0]<<std::endl;
        std::cout<<"FlashGG: Max domain = "<<MinMaxDomain[X][1]<<" "<<MinMaxDomain[Y][1]<<" "<<MinMaxDomain[Z][1]<<std::endl;
        std::cout<<"FlashGG: Length of domain (L) = "<<L[X]<<" "<<L[Y]<<" "<<L[Z]<<std::endl;
        if (grid_type == 'A') {
            std::cout<<"FlashGG: Min cell size = "<<Dmin[X]<<" "<<Dmin[Y]<<" "<<Dmin[Z]<<std::endl;
            std::cout<<"FlashGG: Max cell size = "<<Dmax[X]<<" "<<Dmax[Y]<<" "<<Dmax[Z]<<std::endl;
        }
        if (grid_type == 'U') {
            std::cout<<"FlashGG: Cell size = "<<D[0][X]<<" "<<D[0][Y]<<" "<<D[0][Z]<<std::endl;
        }
        // print by-block info, if keyword set
        if (by_block) {
            for (int b = 0; b < NumBlocks; b++) {
              std::cout<<"FlashGG: cell size (D) = "<<D[b][X]<<" "<<D[b][Y]<<" "<<D[b][Z]<<std::endl;
              std::cout<<"FlashGG: block="<<b<<": Min BBox = "<<BoundingBox[b][X][0]<<" "<<BoundingBox[b][Y][0]<<" "<<BoundingBox[b][Z][0]<<std::endl;
              std::cout<<"FlashGG: block="<<b<<": Max BBox = "<<BoundingBox[b][X][1]<<" "<<BoundingBox[b][Y][1]<<" "<<BoundingBox[b][Z][1]<<std::endl;
              std::cout<<"FlashGG: block="<<b<<": LBlock = "<<LBlock[b][X]<<" "<<LBlock[b][Y]<<" "<<LBlock[b][Z]<<std::endl;
            }
        }
    };

    // Domain decomposition by blocks (overloaded).
    // If AMR: leaf blocks; if UG: pseudo blocks (= normal blocks by default).
    // Inputs: MPI rank (MyPE), total number of MPI ranks (NPE).
    // Return indices of blocks for MyPE.
    public: std::vector<int> GetMyBlocks(const int MyPE, const int NPE)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        std::vector<int> BlockList(0);
        if (grid_type == 'A') {
            BlockList = this->GetLeafBlocks();
        }
        if (grid_type == 'U') {
            for (int i = 0; i < NumBlocks_PB; i++)
                BlockList.push_back(i);
        }
        return this->GetMyBlocks(MyPE, NPE, BlockList);
    }
    // Domain decomposition by blocks.
    // Inputs: MPI rank (MyPE), total number of MPI ranks (NPE), total number of blocks to distribute (nB).
    // Return indices of blocks for MyPE.
    public: std::vector<int> GetMyBlocks(const int MyPE, const int NPE, const std::vector<int> BlockList)
    {
        std::vector<int> MyBlocks(0);
        int DivBlocks = ceil( (double)(BlockList.size()) / (double)(NPE) );
        int NPE_main = BlockList.size() / DivBlocks;
        int ModBlocks = BlockList.size() - NPE_main * DivBlocks;
        if (MyPE < NPE_main) { // (NPE_main) cores get DivBlocks blocks
            for (unsigned int ib = 0; ib < DivBlocks; ib++)
                MyBlocks.push_back(BlockList[MyPE*DivBlocks+ib]);
        }
        if (MyPE == 0) std::cout<<"FlashGG: GetMyBlocks: First "<<NPE_main<<" core(s) carry(ies) "<<DivBlocks<<" block(s) (each)."<<std::endl;
        if ((MyPE == NPE_main) && (ModBlocks > 0)) { // core (NPE_main + 1) gets the rest (ModBlocks)
            for (unsigned int ib = 0; ib < ModBlocks; ib++)
                MyBlocks.push_back(BlockList[NPE_main*DivBlocks+ib]);
            std::cout<<"FlashGG: GetMyBlocks: Core #"<<NPE_main+1<<" carries "<<ModBlocks<<" block(s)."<<std::endl;
        }
        int NPE_in_use = NPE_main; if (ModBlocks > 0) NPE_in_use += 1;
        if ((MyPE == 0) && (NPE_in_use < NPE)) {
            std::cout<<"FlashGG: GetMyBlocks: Warning: non-optimal load balancing; "<<NPE-NPE_in_use<<" core(s) remain(s) idle."<<std::endl;
        }
        if (Debug) {
            std::cout<<" ["<<MyPE<<"] FlashGG: GetMyBlocks: MyBlocks =";
            for (unsigned int ib = 0; ib < MyBlocks.size(); ib++)
                std::cout<<" "<<MyBlocks[ib];
            std::cout<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        return MyBlocks;
    }

    /// GetNumDims
    public: int GetNumDims(void) { return NumDims; };
    /// GetNumBlocks
    public: int GetNumBlocks(void) { return NumBlocks; };
    /// GetNumBlocksRep
    public: int GetNumBlocksRep(void) { return NumBlocksRep; };
    /// GetNumBlocksVector
    public: std::vector<int> GetNumBlocksVector(void) { return NumBlocksIn; };
    /// GetNumCellsInBlock
    public: std::vector<int> GetNumCellsInBlock(void) { return NB; };
    /// GetMinMaxDomain
    public: std::vector< std::vector<double> > GetMinMaxDomain(void) { return MinMaxDomain; };
    /// GetN
    public: std::vector<int> GetN(void) { return N; };
    /// GetL
    public: std::vector<double> GetL(void) { return L; };
    /// GetDblock (return cell size as function of block and dim)
    public: std::vector< std::vector<double> > GetDblock(void) { return D; };
    /// GetD (return cell size as function of dim only; assumes UG)
    public: std::vector<double> GetD(void) {
        if (grid_type == 'U') {
            return D[0];
        } else {
            std::cout<<"FlashGG: Error in GetD(). This is not a uniform grid. Use GetDblock() instead."<<std::endl;
            exit(-1);
        }
    };
    /// GetDmin (return minimum cell size)
    public: std::vector<double> GetDmin(void) { return Dmin; };
    /// GetDmax (return maximum cell size)
    public: std::vector<double> GetDmax(void) { return Dmax; };
    /// GetLBlock
    public: std::vector< std::vector<double> > GetLBlock(void) { return LBlock; };
    /// GetBoundingBox
    public: std::vector< std::vector <std::vector<double> > > GetBoundingBox(void) { return BoundingBox; };
    /// GetNodeType
    public: std::vector<int> GetNodeType(void) { return NodeType; };
    /// Get LeafBlocks (get a list of all active (leaf) blocks)
    public: std::vector<int> GetLeafBlocks(void) {
        std::vector<int> LeafBlocks(0);
        for (int ib=0; ib<NumBlocks; ib++) {
            if (NodeType[ib] == 1) { // LEAF block
                LeafBlocks.push_back(ib);
            }
        }
        return LeafBlocks;
    };
    /// GetNumBlocks_PB
    public: int GetNumBlocks_PB(void) { return NumBlocks_PB; };
    /// GetNumBlocksVector_PB
    public: std::vector<int> GetNumBlocksVector_PB(void) { return NumBlocksIn_PB; };
    /// GetNumCellsInBlock_PB
    public: std::vector<int> GetNumCellsInBlock_PB(void) { return NB_PB; };
    /// GetLBlock_PB
    public: std::vector< std::vector<double> > GetLBlock_PB(void) { return LBlock_PB; };
    /// GetBoundingBox_PB
    public: std::vector< std::vector <std::vector<double> > > GetBoundingBox_PB(void) { return BoundingBox_PB; };

    /// ReadBlockVar (overloaded non-MPI)
    public: FLASH_GG_REAL * ReadBlockVar(const int &block, const std::string datasetname, long int &size)
    {
        size = NB[X]*NB[Y]*NB[Z];
        return this->ReadBlockVar(block, datasetname, MPI_COMM_NULL);
    };
    /// ReadBlockVar (overloaded)
    public: FLASH_GG_REAL * ReadBlockVar(const int &block, const std::string datasetname, long int &size, MPI_Comm comm)
    {
        size = NB[X]*NB[Y]*NB[Z];
        return this->ReadBlockVar(block, datasetname, comm);
    };
    /// ReadBlockVar (overloaded, non-MPI)
    public: FLASH_GG_REAL * ReadBlockVar(const int &block, const std::string datasetname)
    {
        return this->ReadBlockVar(block, datasetname, MPI_COMM_NULL);
    };
    /// ReadBlockVar (overloaded)
    public: FLASH_GG_REAL * ReadBlockVar(const int &block, const std::string datasetname, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r', comm);
        FLASH_GG_REAL * DataPointer = new FLASH_GG_REAL[NB[X]*NB[Y]*NB[Z]];
        hsize_t offset[4] = {block % NumBlocks, 0, 0, 0}; // note that the % NumBlocks takes care of PBCs (if called with a block replica index)
        hsize_t count[4] = {1, NB[Z], NB[Y], NB[X]};
        hsize_t out_offset[3] = {0, 0, 0};
        hsize_t out_count[3] = {NB[Z], NB[Y], NB[X]};
        HDFInput.read_slab(DataPointer, datasetname, FLASH_GG_H5_REAL, offset, count, 3, out_offset, out_count, comm);
        HDFInput.close();
        return DataPointer;
    };

    /// ReadBlockVar_PB (overloaded non-MPI)
    public: FLASH_GG_REAL * ReadBlockVar_PB(const int &block_pb, const std::string datasetname, long int &size)
    {
        size = NB_PB[X]*NB_PB[Y]*NB_PB[Z];
        return this->ReadBlockVar_PB(block_pb, datasetname, MPI_COMM_NULL);
    };
    /// ReadBlockVar_PB (overloaded)
    public: FLASH_GG_REAL * ReadBlockVar_PB(const int &block_pb, const std::string datasetname, long int &size, MPI_Comm comm)
    {
        size = NB_PB[X]*NB_PB[Y]*NB_PB[Z];
        return this->ReadBlockVar_PB(block_pb, datasetname, comm);
    };
    /// ReadBlockVar_PB (overloaded, non-MPI)
    public: FLASH_GG_REAL * ReadBlockVar_PB(const int &block_pb, const std::string datasetname)
    {
        return this->ReadBlockVar_PB(block_pb, datasetname, MPI_COMM_NULL);
    };
    /// ReadBlockVar_PB (overloaded)
    public: FLASH_GG_REAL * ReadBlockVar_PB(const int &block_pb, const std::string datasetname, MPI_Comm comm)
    {
        // create new PB dataset pointer
        FLASH_GG_REAL * DataPointer = new FLASH_GG_REAL[NB_PB[X]*NB_PB[Y]*NB_PB[Z]];
        // open file for read
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r', comm);
        // find all file blocks that overlap this PB
        std::vector<int> AffectedBlocks = this->GetAffectedBlocks(BoundingBox_PB[block_pb]);
        for (unsigned int ib = 0; ib < AffectedBlocks.size(); ib++)
        {
            int block = AffectedBlocks[ib];
            std::vector<int> cb_offset(NumDims); // file block cell offset
            std::vector<int> cb_count(NumDims); // file block cell count
            std::vector<int> cpb_offset(NumDims); // pseudo block cell offset
            std::vector<int> cpb_count(NumDims); // pseudo block cell count
            for (int dim = 0; dim < NumDims; dim++) {
                // if left edge of PB is inside current file block
                if ( (BoundingBox_PB[block_pb][dim][0] > BoundingBox[block][dim][0]) &&
                     (BoundingBox_PB[block_pb][dim][0] < BoundingBox[block][dim][1] ) ) {
                    cb_offset[dim] = (int)((BoundingBox_PB[block_pb][dim][0]+D[block][dim]/2-BoundingBox[block][dim][0])/LBlock[block][dim]*NB[dim]);
                    cpb_offset[dim] = 0;
                } else {
                    cb_offset[dim] = 0;
                    cpb_offset[dim] = (int)((BoundingBox[block][dim][0]+D[block][dim]/2-BoundingBox_PB[block_pb][dim][0])/LBlock_PB[block_pb][dim]*NB_PB[dim]);
                }
                // if right edge of PB is inside current file block
                if ( (BoundingBox_PB[block_pb][dim][1] > BoundingBox[block][dim][0]) &&
                     (BoundingBox_PB[block_pb][dim][1] < BoundingBox[block][dim][1] ) ) {
                    cb_count[dim] = (int)((BoundingBox_PB[block_pb][dim][1]+D[block][dim]/2-BoundingBox[block][dim][0])/LBlock[block][dim]*NB[dim]) - cb_offset[dim];
                } else {
                    cb_count[dim] = NB[dim] - cb_offset[dim];
                }
            }
            // HDF5 offsets and count for slab selections
            hsize_t offset[4] = {block, cb_offset[Z], cb_offset[Y], cb_offset[X]};
            hsize_t count[4] = {1, cb_count[Z], cb_count[Y], cb_count[X]};
            hsize_t out_offset[3] = {cpb_offset[Z], cpb_offset[Y], cpb_offset[X]};
            hsize_t out_count[3] = {cb_count[Z], cb_count[Y], cb_count[X]};
            hsize_t total_out_count[3] = {NB_PB[Z], NB_PB[Y], NB_PB[X]};
            if (Debug) {
                std::cout<<"block = "<<block<<std::endl;
                for (int dim = 0; dim < NumDims; dim++) {
                    std::cout<<"cb_count["<<dim<<"] = "<<cb_count[dim]
                             <<", cb_offset["<<dim<<"] = "<<cb_offset[dim]
                             <<", cpb_offset["<<dim<<"] = "<<cpb_offset[dim]<<std::endl;
                }
            }
            // HDFIO for reading slab
            HDFInput.read_slab(DataPointer, datasetname, FLASH_GG_H5_REAL, offset, count, 3, out_offset, out_count, total_out_count, comm);
        }
        HDFInput.close();
        return DataPointer;
    };

    /// ReadBlockVar_PB_Interpolated
    public: FLASH_GG_REAL * ReadBlockVar_PB_Interpolated(const int &block_pb, const std::string datasetname, const bool mass_weighting)
    {
        // get uniform grid for this PB from file blocks
        return this->GetUniformGrid(NB_PB, BoundingBox_PB[block_pb], datasetname, mass_weighting);
    };

    // GetUniformGrid (overloaded: without PBCs)
    public: FLASH_GG_REAL * GetUniformGrid(const std::vector<int> np,
                                           const std::vector< std::vector<double> > bounds,
                                           const std::string datasetname, const bool mass_weighting)
    {
        return this->GetUniformGrid(np, bounds, datasetname, mass_weighting, false);
    }
    // GetUniformGrid
    public: FLASH_GG_REAL * GetUniformGrid(const std::vector<int> np,
                                           const std::vector< std::vector<double> > bounds,
                                           const std::string datasetname, const bool mass_weighting,
                                           const bool periodic_boundary_conditions)
    {
        if (Debug) std::cout<<"FlashGG: GetUniformGrid: entering."<<std::endl;
        // create GRID3D uniform grid
        assert (np.size() == 3); // this function is currently only implemented for 3D
        assert (bounds.size() == 3); // this function is currently only implemented for 3D
        GRID3D grid_data = GRID3D(np[X], np[Y], np[Z]);
        grid_data.set_bnds(bounds[X][0], bounds[X][1], bounds[Y][0], bounds[Y][1], bounds[Z][0], bounds[Z][1]);
        grid_data.clear();
        // GRID3D for density in case of mass-weighting
        GRID3D grid_dens;
        if (mass_weighting) {
            grid_dens = GRID3D(np[X], np[Y], np[Z]);
            grid_dens.set_bnds(bounds[X][0], bounds[X][1], bounds[Y][0], bounds[Y][1], bounds[Z][0], bounds[Z][1]);
            grid_dens.clear();
        }
        // if PBCs, create block replicas (extend BoundingBox and NumBlocksRep)
        if (periodic_boundary_conditions) this->AddBlockReplicasPBCs();
        // find affected blocks in file
        if (Debug) std::cout<<"FlashGG: GetUniformGrid: finding affected blocks..."<<std::endl;
        // returns all affected block indices (including block replicas if AddBlockReplicasPBCs was called earlier)
        std::vector<int> AffectedBlocks = this->GetAffectedBlocks(bounds);
        // loop over affected blocks
        if (Debug) std::cout<<"FlashGG: GetUniformGrid: looping..."<<std::endl;
        for (unsigned int ib = 0; ib < AffectedBlocks.size(); ib++)
        {
            int b_all = AffectedBlocks[ib];
            int b = b_all % NumBlocks; // take care of PBCs (if present)
            FLASH_GG_REAL * block_data = this->ReadBlockVar(b, datasetname);
            FLASH_GG_REAL * block_dens = 0;
            if (mass_weighting) block_dens = this->ReadBlockVar(b, "dens");
            // loop over cells in this block and assign data to GRID3D
            for (int k = 0; k < NB[Z]; k++)
                for (int j = 0; j < NB[Y]; j++)
                    for (int i = 0; i < NB[X]; i++) {
                        long index = k*NB[X]*NB[Y] + j*NB[X] + i;
                        double dvol = D[b][X]*D[b][Y]*D[b][Z];
                        double data = block_data[index]*dvol;
                        if (mass_weighting) data *= block_dens[index];
                        std::vector<double> cc = this->CellCenter(b_all, i, j, k); /// this function can also take block replica indices
                        grid_data.add_coord_fields(cc[X], cc[Y], cc[Z], D[b][X], D[b][Y], D[b][Z], data);
                        if (mass_weighting) grid_dens.add_coord_fields(cc[X], cc[Y], cc[Z], D[b][X], D[b][Y], D[b][Z], block_dens[index]*dvol);
                    }

            delete [] block_data;
            if (block_dens) delete [] block_dens;
        }
        // if PBCs, remove block replicas (shrink BoundingBox and NumBlocksRep)
        if (periodic_boundary_conditions) this->RemoveBlockReplicasPBCs();
        // prepare output
        if (Debug) std::cout<<"FlashGG: GetUniformGrid: copying to output..."<<std::endl;
        int ntot = grid_data.get_ntot();
        FLASH_GG_REAL * grid_out = new FLASH_GG_REAL[ntot];
        for (int n = 0; n < ntot; n++) {
            if (mass_weighting)
                grid_out[n] = (FLASH_GG_REAL)(grid_data.field[n]/grid_dens.field[n]);
            else
                grid_out[n] = (FLASH_GG_REAL)(grid_data.field[n]);
        }
        if (Debug) std::cout<<"FlashGG: GetUniformGrid: exiting."<<std::endl;
        return grid_out;
    }

    /// GetAffectedBlocks
    public: std::vector<int> GetAffectedBlocks(const std::vector< std::vector<double> > bounds)
    {
        std::vector<int> AffectedBlocks(0);
        assert (bounds.size() == NumDims);
        std::vector<bool> overlap(NumDims);
        for (int b = 0; b < NumBlocksRep; b++) {
            if (NodeType[b % NumBlocks] == 1) { // LEAF block
                // see if there is overlap in any dim
                for (int dim = 0; dim < NumDims; dim++) {
                    if ( ((bounds[dim][0] <= BoundingBox[b][dim][0]) && (BoundingBox[b][dim][0] <  bounds[dim][1])) ||
                         ((bounds[dim][0] <  BoundingBox[b][dim][1]) && (BoundingBox[b][dim][1] <= bounds[dim][1])) ||
                         ((bounds[dim][0] >= BoundingBox[b][dim][0]) && (BoundingBox[b][dim][1] >= bounds[dim][1])) ||
                         ((bounds[dim][0] <  BoundingBox[b][dim][0]) && (BoundingBox[b][dim][1] <  bounds[dim][1]))    ) {
                        overlap[dim] = true;
                    } else {
                        overlap[dim] = false;
                    }
                }
                // check if all elements of overlap are true and if so, append to AffectedBlocks
                if (std::find(overlap.begin(), overlap.end(), false) == overlap.end())
                    AffectedBlocks.push_back(b);
            }
        }
        return AffectedBlocks;
    }

    // AddBlockReplicasPBCs (extend BoundingBox to allow for PBCs)
    public: void AddBlockReplicasPBCs(void)
    {
        // resize (extend) BoundingBox to carry block replicas (always appended after the active blocks, starting at index NumBlocks)
        NumBlocksRep = pow(3, NumDims) * NumBlocks; // reset total number of blocks (now including replicas)
        BoundingBox.resize(NumBlocksRep);
        for (int b_rep = 0; b_rep < NumBlocksRep; b_rep++) {
            BoundingBox[b_rep].resize(NumDims);
            for (int dim = 0; dim < NumDims; dim++) {
                BoundingBox[b_rep][dim].resize(2);
            }
        }
        // generate block replicas for PBCs by setting the BoundingBox'es of the block replicas
        int b_rep_factor = 0;
        int pbc_x_nrep = 1;
        int pbc_y_nrep = 0; if (NumDims > 1) pbc_y_nrep = 1;
        int pbc_z_nrep = 0; if (NumDims > 2) pbc_z_nrep = 1;
        for (int pbc_z = -pbc_z_nrep; pbc_z <= pbc_z_nrep; pbc_z++) { // loop over replicas per dim
            for (int pbc_y = -pbc_y_nrep; pbc_y <= pbc_y_nrep; pbc_y++) { // loop over replicas per dim
                for (int pbc_x = -pbc_x_nrep; pbc_x <= pbc_x_nrep; pbc_x++) { // loop over replicas per dim
                    if ((pbc_x == 0) && (pbc_y == 0) && (pbc_z == 0)) continue; // skip the centre (it's already the original set of blocks)
                    b_rep_factor++;
                    // loop over all active blocks
                    for (int b = 0; b < NumBlocks; b++) {
                        int b_rep = b_rep_factor*NumBlocks + b; // block replica index into BoundingBox (original block index b = b_rep % NumBlocks)
                        assert(b_rep < NumBlocksRep);
                        for (int minmax = 0; minmax < 2; minmax++) {
                            BoundingBox[b_rep][X][minmax] = BoundingBox[b][X][minmax] + pbc_x*L[X];
                            if (NumDims > 1) BoundingBox[b_rep][Y][minmax] = BoundingBox[b][Y][minmax] + pbc_y*L[Y];
                            if (NumDims > 2) BoundingBox[b_rep][Z][minmax] = BoundingBox[b][Z][minmax] + pbc_z*L[Z];
                        } // minmax
                    } // b
                } // pbc_x
            } // pbc_y
        } // pbc_z
    }

    // RemoveBlockReplicasPBCs (shrink BoundingBox to original size)
    public: void RemoveBlockReplicasPBCs(void)
    {
        BoundingBox.resize(NumBlocks);
        NumBlocksRep = NumBlocks;
    }

    /// CreateDataset
    public: void CreateDataset(const std::string datasetname, std::vector<int> Dimensions, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        bool datasetname_exists = false;
        int n_datasets = HDFInput.getNumberOfDatasets();
        for (int n=0; n<n_datasets; n++)
        {
            std::string datasetname_in_file = HDFInput.getDatasetname(n);
            if (datasetname_in_file == datasetname)
            {
                datasetname_exists = true;
            }
        }
        if (datasetname_exists)
        {
            //std::cout<<"FlashGG_mpi: Datasetname '"<<datasetname<<"' already exists in file. SKIPPING creation!"<<std::endl;
        }
        else
        {
            HDFInput.create_dataset(datasetname, Dimensions, FLASH_GG_H5_REAL, comm);
        }
        HDFInput.close();
    };

    /// OverwriteBlockVar (overloaded for typical blocks of this GG)
    public: void OverwriteBlockVar(const int &block, const std::string datasetname, FLASH_GG_REAL * const DataPointer, MPI_Comm comm)
    {
        this->OverwriteBlockVar(block, NB, datasetname, DataPointer, comm);
    };
    /// OverwriteBlockVar (takes the number of cells of that block variable,
    /// in case we created a new block var with different cell dimensions)
    public: void OverwriteBlockVar(const int &block, const std::vector<int> myNB, 
                                   const std::string datasetname, FLASH_GG_REAL * const DataPointer, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        hsize_t offset[4] = {block, 0, 0, 0};
        hsize_t count[4] = {1, myNB[Z], myNB[Y], myNB[X]};
        hsize_t out_offset[3] = {0, 0, 0};
        hsize_t out_count[3] = {myNB[Z], myNB[Y], myNB[X]};
        HDFInput.overwrite_slab(DataPointer, datasetname, FLASH_GG_H5_REAL, offset, count, 3, out_offset, out_count, comm);
        HDFInput.close();
    };

    /// public ReadDatasetNames
    public: std::vector<std::string> ReadDatasetNames(void)
    {
        HDFIO HDFInput = HDFIO();
        std::vector<std::string> ret = HDFInput.getDatasetnames(Inputfilename);
        return ret;
    };

    /// ReadBlockSize
    public: FLASH_GG_REAL * ReadBlockSize(const int &block, long int &size)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        size = NumDims;
        FLASH_GG_REAL * DataPointer = new FLASH_GG_REAL[size];
        hsize_t offset[2] = {block, 0};
        hsize_t count[2] = {1, size};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {size};
        HDFInput.read_slab(DataPointer, "block size", FLASH_GG_H5_REAL, offset, count, 1, out_offset, out_count);
        HDFInput.close();
        return DataPointer;
    };

    /// OverwriteBlockSize
    public: void OverwriteBlockSize(const int &block, FLASH_GG_REAL * const DataPointer, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        hsize_t offset[2] = {block, 0};
        hsize_t count[2] = {1, NumDims};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {NumDims};
        HDFInput.overwrite_slab(DataPointer, "block size", FLASH_GG_H5_REAL, offset, count, 1, out_offset, out_count, comm);
        HDFInput.close();
    };

    /// ReadBoundingBox
    public: FLASH_GG_REAL * ReadBoundingBox(const int &block, long int &size)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        size = NumDims*2;
        FLASH_GG_REAL * DataPointer = new FLASH_GG_REAL[size];
        hsize_t offset[3] = {block, 0, 0};
        hsize_t count[3] = {1, NumDims, 2};
        hsize_t out_offset[2] = {0, 0};
        hsize_t out_count[2] = {NumDims, 2};
        HDFInput.read_slab(DataPointer, bounding_box_datasetname, FLASH_GG_H5_REAL, offset, count, 2, out_offset, out_count);
        HDFInput.close();
        return DataPointer;
    };

    /// OverwriteBoundingBox
    public: void OverwriteBoundingBox(const int &block, FLASH_GG_REAL * const DataPointer, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        hsize_t offset[3] = {block, 0, 0};
        hsize_t count[3] = {1, NumDims, 2};
        hsize_t out_offset[2] = {0, 0};
        hsize_t out_count[2] = {NumDims, 2};
        HDFInput.overwrite_slab(DataPointer, bounding_box_datasetname, FLASH_GG_H5_REAL, offset, count, 2, out_offset, out_count, comm);
        HDFInput.close();
    };

    /// ReadCoordinates
    public: FLASH_GG_REAL * ReadCoordinates(const int &block, long int &size)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        size = NumDims;
        FLASH_GG_REAL * DataPointer = new FLASH_GG_REAL[size];
        hsize_t offset[2] = {block, 0};
        hsize_t count[2] = {1, size};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {size};
        HDFInput.read_slab(DataPointer, "coordinates", FLASH_GG_H5_REAL, offset, count, 1, out_offset, out_count);
        HDFInput.close();
        return DataPointer;
    };

    /// OverwriteCoordinates
    public: void OverwriteCoordinates(const int &block, FLASH_GG_REAL * const DataPointer, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        hsize_t offset[2] = {block, 0};
        hsize_t count[2] = {1, NumDims};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {NumDims};
        HDFInput.overwrite_slab(DataPointer, "coordinates", FLASH_GG_H5_REAL, offset, count, 1, out_offset, out_count, comm);
        HDFInput.close();
    };

    /// ReadGID
    public: int * ReadGID(const int &block, long int &size)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        size = 15;
        int * DataPointer = new int[size];
        hsize_t offset[2] = {block, 0};
        hsize_t count[2] = {1, size};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {size};
        HDFInput.read_slab(DataPointer, "gid", H5T_NATIVE_INT, offset, count, 1, out_offset, out_count);
        HDFInput.close();
        return DataPointer;
    };

    /// OverwriteGID
    public: void OverwriteGID(const int &block, int * const DataPointer, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        hsize_t offset[2] = {block, 0};
        hsize_t count[2] = {1, 15};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {15};
        HDFInput.overwrite_slab(DataPointer, "gid", H5T_NATIVE_INT, offset, count, 1, out_offset, out_count, comm);
        HDFInput.close();
    };

    /// public ReadNodeType
    public: int * ReadNodeType(const int &block, long int &size)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        size = 1;
        int * DataPointer = new int[size];
        hsize_t offset[1] = {block};
        hsize_t count[1] = {size};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {size};
        HDFInput.read_slab(DataPointer, node_type_datasetname, H5T_NATIVE_INT, offset, count, 1, out_offset, out_count);
        HDFInput.close();
        return DataPointer;
    };

    /// OverwriteNodeType
    public: void OverwriteNodeType(const int &block, int * const DataPointer, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        hsize_t offset[1] = {block};
        hsize_t count[1] = {1};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {1};
        HDFInput.overwrite_slab(DataPointer, node_type_datasetname, H5T_NATIVE_INT, offset, count, 1, out_offset, out_count, comm);
        HDFInput.close();
    };

    /// public ReadProcessorNumber
    public: int * ReadProcessorNumber(const int &block, long int &size)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        size = 1;
        int * DataPointer = new int[size];
        hsize_t offset[1] = {block};
        hsize_t count[1] = {size};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {size};
        HDFInput.read_slab(DataPointer, "processor number", H5T_NATIVE_INT, offset, count, 1, out_offset, out_count);
        HDFInput.close();
        return DataPointer;
    };

    /// OverwriteProcessorNumber
    public: void OverwriteProcessorNumber(const int &block, int * const DataPointer, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        hsize_t offset[1] = {block};
        hsize_t count[1] = {1};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {1};
        HDFInput.overwrite_slab(DataPointer, "processor number", H5T_NATIVE_INT, offset, count, 1, out_offset, out_count, comm);
        HDFInput.close();
    };

    /// public ReadRefineLevel
    public: int * ReadRefineLevel(const int &block, long int &size)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        size = 1;
        int * DataPointer = new int[size];
        hsize_t offset[1] = {block};
        hsize_t count[1] = {size};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {size};
        HDFInput.read_slab(DataPointer, "refine level", H5T_NATIVE_INT, offset, count, 1, out_offset, out_count);
        HDFInput.close();
        return DataPointer;
    };

    /// OverwriteRefineLevel
    public: void OverwriteRefineLevel(const int &block, int * const DataPointer, MPI_Comm comm)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'w', comm);
        hsize_t offset[1] = {block};
        hsize_t count[1] = {1};
        hsize_t out_offset[1] = {0};
        hsize_t out_count[1] = {1};
        HDFInput.overwrite_slab(DataPointer, "refine level", H5T_NATIVE_INT, offset, count, 1, out_offset, out_count, comm);
        HDFInput.close();
    };

    /// ReadUnknownNames
    public: std::vector<std::string> ReadUnknownNames(void)
    {
        hid_t File_id = H5Fopen(Inputfilename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t dataset = H5Dopen(File_id, "unknown names", H5P_DEFAULT);
        hid_t dataspace = H5Dget_space(dataset);
        const int rank = H5Sget_simple_extent_ndims(dataspace);
        hsize_t dimens_2d[rank];
        H5Sget_simple_extent_dims(dataspace, dimens_2d, NULL);
        const int n_names = dimens_2d[0];
        // mallocate output
        char ** unk_labels = (char **) malloc (n_names * sizeof (char *));
        // determine whether this is a string of fixed or variable size
        hid_t string_type = H5Tcopy(H5T_C_S1); herr_t status;
        hid_t filetype = H5Dget_type(dataset);
        htri_t variable_length_string = H5Tis_variable_str(filetype);
        // in case it's a string of fixed length:
        if (variable_length_string == 0) {
            size_t string_size = H5Tget_size(filetype); string_size++; /* Make room for null terminator */
            // some horrible additional pointer stuff allocation
            unk_labels[0] = (char *) malloc (n_names * string_size * sizeof (char));
            for (int i = 1; i < n_names; i++) unk_labels[i] = unk_labels[0] + i * string_size;
            status = H5Tset_size(string_type, string_size);
            status = H5Dread(dataset, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, unk_labels[0]);
        }
        // in case it's a string of variable length:
        if (variable_length_string == 1) {
            status = H5Tset_size(string_type, H5T_VARIABLE);
            status = H5Dread(dataset, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, unk_labels);
        }
        // copy c-strings into std::strings
        std::vector<std::string> ret(n_names);
        for (int i = 0; i < n_names; i++) {
            ret[i] = unk_labels[i];
        }
        // garbage collection
        free (unk_labels[0]);
        free (unk_labels);
        H5Tclose(string_type);
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Fclose(File_id);
        return ret;
    };

    /// OverwriteUnknownNames (note that this overwrites STRSIZE 4 with STRSIZE 40)
    /// This must be called by the Master Processor only
    public: void OverwriteUnknownNames(const std::vector<std::string> unknown_names)
    {
        // copy input strings -> c-strings (of fixed size)
        size_t string_size = 40;
        const int n_names = unknown_names.size();
        char unk_labels[n_names][string_size];
        for (int i = 0; i < n_names; i++) std::strcpy(unk_labels[i], unknown_names[i].c_str());
        // open file, delete previous unknown names and write new one
        hid_t File_id = H5Fopen(Inputfilename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        hid_t string_type = H5Tcopy(H5T_C_S1); herr_t status;
        status = H5Tset_size(string_type, string_size);
        H5Ldelete(File_id, "unknown names", H5P_DEFAULT); // delete dataset
        const int rank = 2; hsize_t dimens_2d[rank]; dimens_2d[0] = n_names; dimens_2d[1] = 1;
        hid_t dataspace = H5Screate_simple(rank, dimens_2d, NULL);
        hid_t dataset = H5Dcreate(File_id, "unknown names", string_type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, unk_labels[0]);
        H5Tclose(string_type);
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Fclose(File_id);
    };

    /// BlockIndex (x, y, z)
    public: inline int BlockIndex(const double &x, const double &y, const double &z)
    {
        int ix = (int)((x-MinMaxDomain[X][0])/L[X]*NumBlocksIn[X]);
        int iy = (int)((y-MinMaxDomain[Y][0])/L[Y]*NumBlocksIn[Y]);
        int iz = (int)((z-MinMaxDomain[Z][0])/L[Z]*NumBlocksIn[Z]);
        int block_index = iz*NumBlocksIn[X]*NumBlocksIn[Y] + iy*NumBlocksIn[X] + ix;
        return block_index;
    };
    /// BlockIndex_PB (x, y, z)
    public: inline int BlockIndex_PB(const double &x, const double &y, const double &z)
    {
        int ix = (int)((x-MinMaxDomain[X][0])/L[X]*NumBlocksIn_PB[X]);
        int iy = (int)((y-MinMaxDomain[Y][0])/L[Y]*NumBlocksIn_PB[Y]);
        int iz = (int)((z-MinMaxDomain[Z][0])/L[Z]*NumBlocksIn_PB[Z]);
        int block_index = iz*NumBlocksIn_PB[X]*NumBlocksIn_PB[Y] + iy*NumBlocksIn_PB[X] + ix;
        return block_index;
    };

    /// CellIndexBlock (block, x, y, z)
    public: inline std::vector<int> CellIndexBlock(const int &block, const double &x, const double &y, const double &z)
    {
        std::vector<int> cell_index(3);
        int b = block % NumBlocks; // take care of PBCs if present (here only for LBlock, access orig block ind)
        cell_index[X] = (int)((x-BoundingBox[block][X][0])/LBlock[b][X]*NB[X]);
        cell_index[Y] = (int)((y-BoundingBox[block][Y][0])/LBlock[b][Y]*NB[Y]);
        cell_index[Z] = (int)((z-BoundingBox[block][Z][0])/LBlock[b][Z]*NB[Z]);
        return cell_index;
    };
    /// CellIndexBlock_PB (block, x, y, z)
    public: inline std::vector<int> CellIndexBlock_PB(const int &block, const double &x, const double &y, const double &z)
    {
        std::vector<int> cell_index(3);
        cell_index[X] = (int)((x-BoundingBox_PB[block][X][0])/LBlock_PB[block][X]*NB_PB[X]);
        cell_index[Y] = (int)((y-BoundingBox_PB[block][Y][0])/LBlock_PB[block][Y]*NB_PB[Y]);
        cell_index[Z] = (int)((z-BoundingBox_PB[block][Z][0])/LBlock_PB[block][Z]*NB_PB[Z]);
        return cell_index;
    };

    /// CellIndexDomain (x, y, z)
    public: inline std::vector<int> CellIndexDomain(const double &x, const double &y, const double &z)
    {
        std::vector<int> cell_index(3);
        cell_index[X] = (int)((x-MinMaxDomain[X][0])/L[X]*N[X]);
        cell_index[Y] = (int)((y-MinMaxDomain[Y][0])/L[Y]*N[Y]);
        cell_index[Z] = (int)((z-MinMaxDomain[Z][0])/L[Z]*N[Z]);
        return cell_index;
    };

    /// CellCenter (block index, cell index i, j, k)
    public: inline std::vector<double> CellCenter(const int &block, const int &i, const int &j, const int &k)
    {
        std::vector<double> cell_center(3);
        int b = block % NumBlocks; // take care of PBCs if present (here only for D, access orig block ind)
        cell_center[X] = BoundingBox[block][X][0]+((double)(i)+0.5)*D[b][X];
        cell_center[Y] = BoundingBox[block][Y][0]+((double)(j)+0.5)*D[b][Y];
        cell_center[Z] = BoundingBox[block][Z][0]+((double)(k)+0.5)*D[b][Z];
        return cell_center;
    };
    /// CellCenter_PB (block index, cell index i, j, k)
    public: inline std::vector<double> CellCenter_PB(const int &block, const int &i, const int &j, const int &k)
    {
        std::vector<double> cell_center(3);
        cell_center[X] = BoundingBox_PB[block][X][0]+((double)(i)+0.5)*D[0][X];
        cell_center[Y] = BoundingBox_PB[block][Y][0]+((double)(j)+0.5)*D[0][Y];
        cell_center[Z] = BoundingBox_PB[block][Z][0]+((double)(k)+0.5)*D[0][Z];
        return cell_center;
    };

    /// CellCenter (block index, cell index)
    public: inline std::vector<double> CellCenter(const int &block, const long &cellindex)
    {
        std::vector<double> cell_center(3);
        int kmod = cellindex % NBXY;
        int k = cellindex / NBXY;
        int j = kmod / NB[X];
        int i = kmod % NB[X];
        int b = block % NumBlocks; // take care of PBCs if present (here only for D, access orig block ind)
        cell_center[X] = BoundingBox[block][X][0]+((double)(i)+0.5)*D[b][X];
        cell_center[Y] = BoundingBox[block][Y][0]+((double)(j)+0.5)*D[b][Y];
        cell_center[Z] = BoundingBox[block][Z][0]+((double)(k)+0.5)*D[b][Z];
        return cell_center;
    };
    /// CellCenter_PB (block index, cell index)
    public: inline std::vector<double> CellCenter_PB(const int &block, const long &cellindex)
    {
        std::vector<double> cell_center(3);
        int kmod = cellindex % NBXY_PB;
        int k = cellindex / NBXY_PB;
        int j = kmod / NB_PB[X];
        int i = kmod % NB_PB[X];
        cell_center[X] = BoundingBox_PB[block][X][0]+((double)(i)+0.5)*D[0][X];
        cell_center[Y] = BoundingBox_PB[block][Y][0]+((double)(j)+0.5)*D[0][Y];
        cell_center[Z] = BoundingBox_PB[block][Z][0]+((double)(k)+0.5)*D[0][Z];
        return cell_center;
    };

    /// ReadNumBlocks
    private: void ReadNumBlocks(void)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        std::vector<int>Dim(3);
        Dim = HDFInput.getDims(bounding_box_datasetname);
        HDFInput.close();
        NumBlocks    = Dim[0];
        NumBlocksRep = NumBlocks; // default is no block replicas
        NumDims      = Dim[1];
        assert(Dim[2] == 2); // min, max
    };

    /// ReadNumCellsInBlock
    private: void ReadNumCellsInBlock(void)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        std::vector<int>Dim(4);
        Dim = HDFInput.getDims("dens");
        HDFInput.close();
        NB.resize(NumDims);
        NB[X] = Dim[3];
        NB[Y] = Dim[2];
        NB[Z] = Dim[1];
        NBXY = NB[X]*NB[Y];
    };

    /// ReadBoundingBoxAndMinMaxDomain
    private: void ReadBoundingBoxAndMinMaxDomain(void)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        std::vector<int>Dim(3);
        Dim = HDFInput.getDims(bounding_box_datasetname);
        NumBlocks = Dim[0];
        NumDims   = Dim[1];
        assert(Dim[2] == 2); // min, max
        FLASH_GG_REAL * BoundingBoxPointer = new FLASH_GG_REAL[NumBlocks*NumDims*2];
        if (Debug) std::cout << "FlashGG: ReadMinMaxDomain: " << bounding_box_datasetname << std::endl;
        HDFInput.read(BoundingBoxPointer, bounding_box_datasetname, FLASH_GG_H5_REAL);
        HDFInput.close();
        if (Debug) {
            for (int ind = 0; ind < NumBlocks*NumDims*2; ind++) {
                std::cout << "FlashGG:  BoundingBoxPointer["<<ind<<"] = "<<BoundingBoxPointer[ind] << std::endl;
            }
        }
        MinMaxDomain.resize(NumDims);
        for (int dim = 0; dim < NumDims; dim++) {
          MinMaxDomain[dim].resize(2);
          MinMaxDomain[dim][0] = BoundingBoxPointer[2*dim+0]; //init
          MinMaxDomain[dim][1] = BoundingBoxPointer[2*dim+1]; //init
        }
        BoundingBox.resize(NumBlocks);
        LBlock.resize(NumBlocks);
        D.resize(NumBlocks);
        Dmin.resize(NumDims); Dmin[X] = +1e99; Dmin[Y] = +1e99; Dmin[Z] = +1e99;
        Dmax.resize(NumDims); Dmax[X] = -1e99; Dmax[Y] = -1e99; Dmax[Z] = -1e99;
        if (Debug) std::cout << "FlashGG: NumBlocks: " << NumBlocks << std::endl;
        for (int block = 0; block < NumBlocks; block++) {
          BoundingBox[block].resize(NumDims);
          LBlock[block].resize(NumDims);
          D[block].resize(NumDims);
          for (int dim = 0; dim < NumDims; dim++) {
            BoundingBox[block][dim].resize(2);
            for (int minmax = 0; minmax < 2; minmax++) {
                int index = NumDims*2*block + 2*dim + minmax;
                if (Debug) {
                    std::cout << "FlashGG: BoundingBoxPointer["<<index<<
                        "] (block="<<block<<" dim="<<dim<<" minmax="<<minmax<<") = "<<
                        BoundingBoxPointer[index] << std::endl;
                }
                BoundingBox[block][dim][minmax] = BoundingBoxPointer[index];
                if (BoundingBox[block][dim][minmax] < MinMaxDomain[dim][0])
                  MinMaxDomain[dim][0] = BoundingBox[block][dim][minmax];
                if (BoundingBox[block][dim][minmax] > MinMaxDomain[dim][1])
                  MinMaxDomain[dim][1] = BoundingBox[block][dim][minmax];
            }
            LBlock[block][dim] = BoundingBox[block][dim][1]-BoundingBox[block][dim][0];
            D[block][dim] = LBlock[block][dim]/(double)(NB[dim]);
            if (D[block][dim] < Dmin[dim]) Dmin[dim] = D[block][dim];
            if (D[block][dim] > Dmax[dim]) Dmax[dim] = D[block][dim];
          }
        }
        // Check whether we are dealing with a 1D or 2D simulation,
        // in which case the bounding box min and max are the same.
        // In such cases, we reset the bounding box to a finte size,
        // but the number of cells remains 1 in the respective dimension(s)
        for (int block = 0; block < NumBlocks; block++) {
          for (int dim = 0; dim < NumDims; dim++) {
            if (BoundingBox[block][dim][0]==BoundingBox[block][dim][1])
            {
                BoundingBox[block][dim][0] = 0.0;
                BoundingBox[block][dim][1] = 1.0;
                MinMaxDomain[dim][0] = 0.0;
                MinMaxDomain[dim][1] = 1.0;
                LBlock[block][dim] = 1.0;
            }
          }
        }
        L.resize(NumDims);
        NumBlocksIn.resize(NumDims);
        N.resize(NumDims);
        for (int dim = 0; dim < NumDims; dim++) {
          L[dim] = MinMaxDomain[dim][1]-MinMaxDomain[dim][0];
          NumBlocksIn[dim] = (int)(L[dim]/LBlock[0][dim]+0.1); // blocks have same size in UG
          N[dim] = (int)(L[dim]/Dmin[dim]+0.1); // effective maximum resolution
        }
        delete [] BoundingBoxPointer;
    };

    /// private ReadNodeType
    private: void ReadNodeType(void)
    {
        HDFIO HDFInput = HDFIO();
        HDFInput.open(Inputfilename, 'r');
        std::vector<int>Dim(1);
        Dim = HDFInput.getDims(node_type_datasetname);
        NumBlocks = Dim[0];
        int * NodeTypePointer = new int[NumBlocks];
        HDFInput.read(NodeTypePointer, node_type_datasetname, H5T_NATIVE_INT);
        HDFInput.close();
        NodeType.resize(NumBlocks);
        for (int block = 0; block < NumBlocks; block++)
            NodeType[block] = NodeTypePointer[block];
        delete [] NodeTypePointer;
    };

    /// GetGridType
    public: char GetGridType(void)
    {
        std::map<std::string, int> integer_params = this->ReadIntegerParameters();
        char grid_type = 'U';
        if (integer_params.count("lrefine_max") == 1) { grid_type = 'A'; }
        return grid_type;
    };

    /// ReadIntegerParameters
    public: std::map<std::string, int> ReadIntegerParameters(void)
    {
        QUICKFLASH_HDF5 qfh5 = QUICKFLASH_HDF5();
        std::map<std::string, int> integer_params;
        qfh5.read_integer_parameters(Inputfilename, integer_params);
        return integer_params;
    };
    /// ReadRealParameters
    public: std::map<std::string, double> ReadRealParameters(void)
    {
        QUICKFLASH_HDF5 qfh5 = QUICKFLASH_HDF5();
        std::map<std::string, double> real_params;
        qfh5.read_real_parameters(Inputfilename, real_params);
        return real_params;
    };
    /// ReadStringParameters
    public: std::map<std::string, std::string> ReadStringParameters(void)
    {
        QUICKFLASH_HDF5 qfh5 = QUICKFLASH_HDF5();
        std::map<std::string, std::string> string_params;
        qfh5.read_string_parameters(Inputfilename, string_params);
        return string_params;
    };
    /// ReadLogicalParameters
    public: std::map<std::string, bool> ReadLogicalParameters(void)
    {
        QUICKFLASH_HDF5 qfh5 = QUICKFLASH_HDF5();
        std::map<std::string, bool> logical_params;
        qfh5.read_logical_parameters(Inputfilename, logical_params);
        return logical_params;
    };

    /// ReadIntegerScalars
    public: std::map<std::string, int> ReadIntegerScalars(void)
    {
        QUICKFLASH_HDF5 qfh5 = QUICKFLASH_HDF5();
        std::map<std::string, int> integer_scalars;
        qfh5.read_integer_scalars(Inputfilename, integer_scalars);
        return integer_scalars;
    };
    /// ReadRealScalars
    public: std::map<std::string, double> ReadRealScalars(void)
    {
        QUICKFLASH_HDF5 qfh5 = QUICKFLASH_HDF5();
        std::map<std::string, double> real_scalars;
        qfh5.read_real_scalars(Inputfilename, real_scalars);
        return real_scalars;
    };
    /// ReadStringScalars
    public: std::map<std::string, std::string> ReadStringScalars(void)
    {
        QUICKFLASH_HDF5 qfh5 = QUICKFLASH_HDF5();
        std::map<std::string, std::string> string_scalars;
        qfh5.read_string_scalars(Inputfilename, string_scalars);
        return string_scalars;
    };
    /// ReadLogicalScalars
    public: std::map<std::string, bool> ReadLogicalScalars(void)
    {
        QUICKFLASH_HDF5 qfh5 = QUICKFLASH_HDF5();
        std::map<std::string, bool> logical_scalars;
        qfh5.read_logical_scalars(Inputfilename, logical_scalars);
        return logical_scalars;
    };

    public: void OverwriteIntegerScalarsAndRuntimeParameters(void)
    {
        int MyPE = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &MyPE);
        if (MyPE != 0) return;
        std::vector<std::string> fieldnames; std::vector<int> fieldvalues;
        fieldnames.push_back("iprocs"); fieldvalues.push_back(NumBlocksIn[X]);
        fieldnames.push_back("jprocs"); fieldvalues.push_back(NumBlocksIn[Y]);
        fieldnames.push_back("kprocs"); fieldvalues.push_back(NumBlocksIn[Z]);
        this->OverwriteMultipleIntegerScalarsOrRuntimeParameters("integer runtime parameters", fieldnames, fieldvalues);
        fieldnames.push_back("nxb"); fieldvalues.push_back(NB[X]);
        fieldnames.push_back("nyb"); fieldvalues.push_back(NB[Y]);
        fieldnames.push_back("nzb"); fieldvalues.push_back(NB[Z]);
        fieldnames.push_back("globalnumblocks"); fieldvalues.push_back(NumBlocks);
        fieldnames.push_back("splitnumblocks"); fieldvalues.push_back(NumBlocks);
        this->OverwriteMultipleIntegerScalarsOrRuntimeParameters("integer scalars", fieldnames, fieldvalues);
    };

    public: void OverwriteMultipleIntegerScalarsOrRuntimeParameters(std::string datasetname,
                                                                    std::vector<std::string> fieldnames,
                                                                    std::vector<int> fieldvalues)
    {
        for (int i = 0; i < fieldnames.size(); i++) {
            this->OverwriteIntegerScalarOrRuntimeParameter(datasetname, fieldnames[i], fieldvalues[i]);
        }
    };

    public: void OverwriteIntegerScalarOrRuntimeParameter(std::string datasetname, std::string fieldname, int fieldvalue)
    {
        const bool Debug = false;
        if ((datasetname != "integer scalars") && (datasetname != "integer runtime parameters")) {
            std::cout<<"ERROR in call to FlashGG: OverwriteIntegerScalarOrRuntimeParameter"<<std::endl;
        }
        hid_t file_id = H5Fopen(Inputfilename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        hid_t intScalarsId = H5Dopen(file_id, datasetname.c_str(), H5P_DEFAULT);
        hid_t spaceId = H5Dget_space(intScalarsId);
        hsize_t scalarDims[1];
        hsize_t ndims = H5Sget_simple_extent_dims(spaceId, scalarDims, NULL); (void) ndims;
        int nScalars = scalarDims[0];
        struct IntegerScalars{ char name[20]; int value; };
        hid_t datatype = H5Tcreate(H5T_COMPOUND, sizeof(IntegerScalars));
        hid_t string20 = H5Tcopy(H5T_C_S1);
        H5Tset_size(string20, 20);
        H5Tinsert(datatype, "name", HOFFSET(IntegerScalars,name), string20);
        H5Tinsert(datatype, "value", HOFFSET(IntegerScalars,value), H5T_NATIVE_INT);
        IntegerScalars *is = new IntegerScalars[nScalars];
        H5Dread(intScalarsId, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, is);
        for (int i = 0; i < nScalars; i++)
        { 
            if (strncmp(is[i].name, fieldname.c_str(), fieldname.size()) == 0) {
                if (Debug) std::cout<<"FlashGG: Overwriting "<<fieldname<<" of "<<is[i].value<<" with "<<
                    fieldvalue<<" in '"<<datasetname<<"'"<<std::endl;
                is[i].value = fieldvalue;
            }
        }
        H5Dwrite(intScalarsId, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, is);
        H5Tclose(string20);
        H5Tclose(datatype);
        H5Sclose(spaceId);
        H5Dclose(intScalarsId);
        H5Fclose(file_id);
        delete [] is;
    };

}; // end: FlashGG
#endif
