from hhfunctions import HHSEARCH_Functions

if __name__ == '__main__':
    #x = HHSEARCH_Functions('/home/bgovi/Workspace/SPRING/spring/spring/spring.test-data/NP_000028.3-YP_009724389.1/NP_000028.3.hhr' )
    x = HHSEARCH_Functions('/home/bgovi/Workspace/SPRING/spring/spring/spring.test-data/NP_000028.3-YP_009724389.1/NP_000028.3.hhr' )
    #x = HHSEARCH_Functions('/home/bgovi/Workspace/SPRING/spring/spring/spring.test-data/NP_000028.3-YP_009724389.1/12asA.hhr' )
    x.MakeModel('/home/bgovi/Workspace/SPRING/spring/spring/spring.test-data/springdb/chains/1U/1UOH_A.pdb','6M3P_E')
    x.MakeModel('/home/bgovi/Workspace/SPRING/spring/spring/spring.test-data/springdb/chains/6M/6M3P_E.pdb','6M3P_E')


    x = HHSEARCH_Functions('/home/bgovi/Workspace/SPRING/spring/spring/spring.test-data/NP_000028.3-YP_009724389.1/YP_009724389.1.hhr' )
    x.MakeModel('/home/bgovi/Workspace/SPRING/spring/spring/spring.test-data/springdb/chains/7C/7C2K_A.pdb','7C2K_A')
