def echo_input(jt, sup, matl, sect, elem, load): 
    ei = "---USER INPUT---\n\n"

    ei += "JOINT      X-COORD        Y-COORD"
    for i in range(len(jt)):
        ei += "\n{0:^5d}     {1:^-10.2f}     {2:^-10.2f}".format(i+1,jt[i][0],jt[i][1])


    ei += "\n\nSUPPORT    JOINT    U1    U2"
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

def echo_output(gdisp, mlefl, reac):
    pass
