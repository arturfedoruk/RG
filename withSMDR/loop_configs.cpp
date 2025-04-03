// loop configurations for the functions below: 
// {Mt, Mt, Mh, MZ, MW, GFermi, QCDQED_at_MZ, mbmb} 
// {y, ?, z, x, x, z, x, y} -1 loop
// (Mt has two loop configurations)
float config_111111[9] = {0, 0, 0, 0, 0, 0, 1, 1} ; // for QCDQED_at_MZ & mbmb loop 0 doesn't exist
float config_222222[9] = {1, 1, 1, 1, 1, 1, 1, 1} ;
float config_333333[9] = {2, 2, 2, 2, 2, 2, 2, 2} ;
float config_333221[9] = {2, 2, 0, 2, 2, 0, 2, 1} ;
float config_444332[9] = {3, 2, 1, 2.5, 2.5, 1, 3, 2} ; 
float config_444333[9] = {3, 2, 2, 2.5, 2.5, 2, 3, 2} ; // hfor MZ & MW loop 3 doesn't exist

const int nconfigs = 6;
float* loop_configs[nconfigs] = {config_111111, config_222222, config_333333, config_333221, config_444332, config_444333};
string loop_names[nconfigs] = {"(1,1,1;1,1;1)","(2,2,2;2,2;2)","(3,3,3;3,3;3)","(3,3,3;2,2;1)","(4,4,4;3,3;2)","(4,4,4;3,3;3)"};
