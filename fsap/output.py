def echo_input(jt, sup, matl, sect, elem, load): 
    ei = "---USER INPUT---\n\n"

    ei += "JOINT      X-COORD        Y-COORD"
    for i in range(len(jt)):
        ei += "\n{0:^5d}     {1:^-10.2f}     {2:^-10.2f}".format(i+1,jt[i][0],jt[i][1])


    ei += "\n\nSUPPORT    JOINT    UX    UY"
    for i in range(len(sup)):
        ei += "\n{0:^7d}    {1:^5d}    {2:^2d}    {3:^2d}".format(i+1,
                                                                  sup[i][0],
                                                                  sup[i][1],
                                                                  sup[i][2])

    ei += "\n\nMATERIAL    E"
    for i in range(len(matl)):
        ei += "\n{0:^8d}    {1:^1.1f}".format(i+1, matl[i])

    ei += "\n\nSECTION    A"
    for i in range(len(sect)):
        ei += "\n{0:^7d}    {1:^1.1f}".format(i+1, sect[i])

    ei += "\n\nMEMBER    BEGJT    ENDJT    MATL    SECT"
    for i in range(len(elem)):
        ei += "\n{0:^6d}    {1:^5d}    {2:^5d}    {3:^4d}    {4:^4d}".format(i+1,
                                                             elem[i][0],
                                                             elem[i][1],
                                                             elem[i][2],
                                                             elem[i][3])

    ei += "\n\nLOAD    JOINT    FX       FY"
    for i in range(len(load)):
        ei += "\n{0:^4d}    {1:^5d}    {2:^2.3f}    {3:^2.3f}".format(i+1,
                                                                      load[i][0],
                                                                      load[i][1],
                                                                      load[i][2])

    print ei

def echo_output(gdisp, mlefl, reac, scv, jt, num_scj, num_dof):
    eo = "\n\n---PROGRAM OUTPUT---\n\n"

    eo += "JOINT DISPLACEMENTS"
    eo += "\n{0:^5s}      {1:^6s}      {2:^6s}".format("JOINT","UX","UY")
    for jt_num in range(1,len(jt)+1):
        str_coord_index_beg = (jt_num - 1)*num_scj
        str_coord_beg = scv[str_coord_index_beg ]
        str_coord_index_end = (jt_num - 1)*num_scj + 1
        str_coord_end = scv[str_coord_index_end ]
        if str_coord_beg <= num_dof:
            u1 = gdisp[str_coord_beg - 1][0]
        else:
            u1 = 0.0

        if str_coord_end <= num_dof:
            u2 = gdisp[str_coord_end - 1][0]
        else:
            u2 = 0.0

        eo += "\n{0:^5d}      {1:^ 6.3f}      {2:^ 6.3f}".format(jt_num, u1, u2)


    eo += "\n\nMEMBER FORCES"
    eo += "\n{0:^6s}    {1:^10s}".format("MEMBER","Fx")
    for m in range(1,len(mlefl)+1):
        eo += "\n{0:^6d}    {1:^ 10.3f}".format(m,mlefl[m-1][2])

    
    eo += "\n\nSUPPORT REACTIONS"
    eo += "\n{0:^5s}      {1:^10s}      {2:^10s}".format("JOINT","FX","FY")
    for jt_num in range(1, len(jt)+1):
        str_coord_index_beg = (jt_num - 1)*num_scj
        str_coord_beg = scv[str_coord_index_beg ]
        str_coord_index_end = (jt_num - 1)*num_scj + 1
        str_coord_end = scv[str_coord_index_end ]
        if str_coord_beg > num_dof:
            fx = reac[str_coord_beg - 1 - num_dof]
        else:
            fx = 0.0

        if str_coord_end > num_dof:
            fy = reac[str_coord_end - 1 - num_dof]
        else:
            fy = 0.0

        eo += "\n{0:^5d}      {1:^ 10.3f}      {2:^ 10.3f}".format(jt_num, fx, fy)

    print eo
