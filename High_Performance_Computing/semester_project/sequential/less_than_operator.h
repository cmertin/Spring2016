;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.

    inline bool TreeNode::operator<(TreeNode const &other) const {
#ifdef __DEBUG_TN__
        if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
      std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
      std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
      assert(false);
    }
#endif

       if ((this->m_uiX == other.m_uiX) && (this->m_uiY == other.m_uiY) && (this->m_uiZ == other.m_uiZ)) {
            return ((this->m_uiLevel & ot::TreeNode::MAX_LEVEL) < (other.m_uiLevel & ot::TreeNode::MAX_LEVEL));
        } //end if

#ifdef HILBERT_ORDERING // To use hilbert, put "#define HILBERT_ORDERING" in the code
#ifdef USE_NCA_PROPERTY
        // #pragma message "Hilbert NCA"
        // If you need to initilize the Hilbert table and the rotations for 2D you need to define DENDRO_DIM2
        // Default initialization for 3D case.
        // NOTE: To work the Hilbert Ordering You need the Hilbert Table Initialized.

        unsigned int x1 = m_uiX;
        unsigned int x2 = other.getX();

        unsigned int y1 = m_uiY;
        unsigned int y2 = other.getY();

        unsigned int z1 = m_uiZ;
        unsigned int z2 = other.getZ();

        unsigned len;
        unsigned int maxDepth = m_uiMaxDepth;

        if(this->getLevel()>other.getLevel())
        {
            len=1u<<(this->getMaxDepth()-other.getLevel());
            if(!((x1<x2 || x1>=(x2+len)) || (y1<y2 || y1>=(y2+len)) ||(z1<z2 || z1>=(z2+len))))
                return false;


        }else if(this->getLevel()<other.getLevel())
        {
            len=1u<<(this->getMaxDepth()-this->getLevel());
            if(!((x2<x1 || x2>=(x1+len))||(y2<y1 || y2>=(y1+len))||(z2<z1 || z2>=(z1+len))))
                return true;


        }

        unsigned int maxDiff = (unsigned int)(std::max((std::max((x1^x2),(y1^y2))),(z1^z2)));
        int dim=m_uiDim;

        unsigned int maxDiffBinLen = binOp::binLength(maxDiff);
        //Eliminate the last maxDiffBinLen bits.
        unsigned int ncaX = ((x1>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaY = ((y1>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaZ = ((z1>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaLev = (maxDepth - maxDiffBinLen);

//        if(ncaLev>std::min(this->getLevel(),other.getLevel()))
//        {
//
//            std::cout<<"P1:"<<*(this)<<"\t P2:"<<other<<std::endl;
//            std::cout<<"MaxDiff:"<<maxDiff<<std::endl;
//            std::cout<<"MaxDiffLen:"<<maxDiffBinLen<<std::endl;
//            std::cout<<"NCALEV:"<<ncaLev<<"\t this_lev:"<<this->getLevel()<<"\t other_lev:"<<other.getLevel()<<std::endl;
//
//            std::cout<<std::endl;
//
//        }

        unsigned int index1=0;
        unsigned int index2=0;
        unsigned int num_children=1u<<dim; // This is basically the hilbert table offset
        unsigned int rot_offset=num_children<<1;
        char index_temp=0;
        int current_rot=0;

        //unsigned int b_x,b_y,b_z;
        //unsigned int a,b,c;
        unsigned int mid_bit=m_uiMaxDepth;

        for(int i=0; i<ncaLev;i++)
        {
            mid_bit=m_uiMaxDepth-i-1;

            //b_x=((ncaX&(1<<mid_bit))>>mid_bit);
            //b_y=((ncaY&(1<<mid_bit))>>mid_bit);
            //b_z=((ncaZ&(1<<mid_bit))>>mid_bit);

            // index1=(b_z<<2) + ((b_x^b_z)<<1) + (b_x^b_y^b_z);
            index1= (((ncaZ&(1<<mid_bit))>>mid_bit)<<2)|( (((ncaX&(1<<mid_bit))>>mid_bit)^((ncaZ&(1<<mid_bit))>>mid_bit)) <<1)|(((ncaX&(1<<mid_bit))>>mid_bit)^((ncaY&(1<<mid_bit))>>mid_bit)^((ncaZ&(1<<mid_bit))>>mid_bit));
            //index_temp=rotations[rot_offset*current_rot+num_children+index1]-'0';
            current_rot=HILBERT_TABLE[current_rot*num_children+index1];

        }


        mid_bit--;
        index1= (((z1&(1<<mid_bit))>>mid_bit)<<2)|( (((x1&(1<<mid_bit))>>mid_bit)^((z1&(1<<mid_bit))>>mid_bit)) <<1)|(((x1&(1<<mid_bit))>>mid_bit)^((y1&(1<<mid_bit))>>mid_bit)^((z1&(1<<mid_bit))>>mid_bit));
        index2= (((z2&(1<<mid_bit))>>mid_bit)<<2)|( (((x2&(1<<mid_bit))>>mid_bit)^((z2&(1<<mid_bit))>>mid_bit)) <<1)|(((x2&(1<<mid_bit))>>mid_bit)^((y2&(1<<mid_bit))>>mid_bit)^((z2&(1<<mid_bit))>>mid_bit));


        return rotations[rot_offset*current_rot+num_children+index1] < rotations[rot_offset*current_rot+num_children+index2];




#else
        //#pragma message "Hilbert"
        //NOTE: We can remove this code later. // WARNING: External Variable set is mandatory.
        G_MAX_DEPTH=m_uiMaxDepth;
        G_dim=m_uiDim;
	    return hilbert_order(p1,p2);
#endif
#else
        //#ifdef USE_NCA_PROPERTY

	      //  #pragma message "Morton NCA"
	      //  return morton_order_NCA(p1,p2);
        // #else
	    //#pragma message "Morton"
          // -- original Morton
          // first compare the x, y, and z to determine which one dominates ...
          //Ancestor is smaller.
          if ((this->m_uiX == other.m_uiX) && (this->m_uiY == other.m_uiY) && (this->m_uiZ == other.m_uiZ)) {
            return ((this->m_uiLevel & ot::TreeNode::MAX_LEVEL) < (other.m_uiLevel & ot::TreeNode::MAX_LEVEL));
          } //end if

          unsigned int x = (m_uiX ^ other.m_uiX);
          unsigned int y = (m_uiY ^ other.m_uiY);
          unsigned int z = (m_uiZ ^ other.m_uiZ);

          //Default pref: z > y > x.
          unsigned int maxC = z;
          unsigned int yOrx = y;
          if (yOrx < x) {if ((x ^ yOrx) >= yOrx) {yOrx = x;}
          }
          if (maxC < yOrx) {if ((maxC ^ yOrx) >= maxC) {maxC = yOrx;}
          }

          if (maxC == z) {return (m_uiZ < other.m_uiZ); } else if (maxC == y) {return (m_uiY < other.m_uiY); } else {return (m_uiX < other.m_uiX); }
        // -- original Morton

    // #endif
#endif

    } //end function
