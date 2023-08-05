

import matplotlib.pyplot as plt


def visualize_teleportation_cost_k_equals_2():

    circuits = ['Sym9_147','4gt5_76','4mod7','Figure4','Parity_247','QFT4','QFT8','QFT16']
    zomorodi_tc = [48,14,10,4,2,8,38,133]
    #kl_tc = [14,6,4,2,2,4,12,32]
    #sp_tc = [18,10,4,2,2,4,8,26]
    ga_tc = [12,4,4,2,2,4,8,16]
    ga_kl_tc = [6,4,4,2,2,4,8,60]
    ga_sp_tc = [12,4,4,2,2,4,8,60]

    x_indexes = range(len(circuits))
    width = 0.2
    plt.bar(x_indexes,zomorodi_tc,width=width,label='TK_GA_Zomorodi')
    plt.bar([i + width for i in x_indexes],ga_kl_tc,width=width,label='TK_GA_KL')
    plt.bar([i + width * 2 for i in x_indexes],ga_sp_tc,width=width,label='TK_GA_SP')
    plt.bar([i + width * 3 for i in x_indexes],ga_tc,width=width,label='TK_GA')

    for i,v in enumerate(zomorodi_tc):
        plt.text(i - width / 6,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')
    for i,v in enumerate(ga_kl_tc):
        plt.text(i + width ,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')
    for i,v in enumerate(ga_sp_tc):
        plt.text(i + width * 2.2,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')
    for i,v in enumerate(ga_tc):
        plt.text(i + width * 3.5,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')

    plt.ylabel('Teleportationskosten')
    plt.title('Vergleich der Teleportationskosten für K=2')
    plt.xticks([i + width  for i in x_indexes],circuits)
    plt.xticks(rotation=45)
    plt.legend()
    plt.show()

def visualize_teleportation_cost_k_equals_3():

    circuits = ['Sym9_147','Rd73_140','4gt5_76','4mod7','Figure4','Parity_247','QFT4','QFT8','QFT16']
    kl_tc = [26,36,30,42,10,14,6,30,108]
    sp_tc = [50,46,20,24,10,20,10,32,80]
    ga_tc = [34,18,10,10,4,4,6,14,30]
    ga_kl_tc = [20,20,10,10,6,4,6,24,110]
    ga_sp_tc = [26,18,6,2,6,6,6,24,120]

    fig,axs = plt.subplots(nrows=1,ncols=5,figsize=(30,10))
    fig.suptitle('Vergleich der Teleportationskosten für K=3')

    width = 0.7
    axs[0].bar(range(len(circuits)),kl_tc,width=width)
    axs[0].set_title('TK_KL')
    for i,v in enumerate(kl_tc):
        axs[0].text(i,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')
    axs[0].set_xticks(range(len(circuits)))
    axs[0].set_xticklabels(circuits,rotation=65,ha='right',fontsize='small')

    axs[1].bar(range(len(circuits)),sp_tc,width=width)
    axs[1].set_title('TK_SP')
    for i,v in enumerate(sp_tc):
        axs[1].text(i,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')
    axs[1].set_xticks(range(len(circuits)))
    axs[1].set_xticklabels(circuits,rotation=65,ha='right',fontsize='small')

    axs[2].bar(range(len(circuits)),ga_kl_tc,width=width)
    axs[2].set_title('TK_GA_KL')
    for i,v in enumerate(ga_kl_tc):
        axs[2].text(i,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')
    axs[2].set_xticks(range(len(circuits)))
    axs[2].set_xticklabels(circuits,rotation=65,ha='right',fontsize='small')

    axs[3].bar(range(len(circuits)),ga_sp_tc,width=width)
    axs[3].set_title('TK_GA_SP')
    for i,v in enumerate(ga_sp_tc):
        axs[3].text(i,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')
    axs[3].set_xticks(range(len(circuits)))
    axs[3].set_xticklabels(circuits,rotation=65,ha='right',fontsize='small')

    axs[4].bar(range(len(circuits)),ga_tc,width=width)
    axs[4].set_title('TK_GA')
    for i,v in enumerate(ga_tc):
        axs[4].text(i,v + 1,str(v),color='black',fontweight='bold',fontsize='small',ha='center')
    axs[4].set_xticks(range(len(circuits)))
    axs[4].set_xticklabels(circuits,rotation=65,ha='right',fontsize='small')

    plt.show()

def visualize_teleportation_cost_k_equals_4():
    circuits = ['Sym9_147','Rd73_140','4gt5_76','4mod7','Figure4','Parity_247','QFT4','QFT8','QFT16']
    kl_tc =    [60,36,30,50,12,16,12,48,158]
    sp_tc =    [68,64,22,50,12,20,12,48,124]
    ga_tc =    [52,32,22,32,12,6,12,24,52]
    ga_kl_tc = [48,20,10,14,8,6,12,30,162]
    ga_sp_tc = [48,34,14,14,8,6,12,34,162]

    fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(30,10))
    fig.suptitle('Vergleich der Teleportationskosten für K=4')

    width = 0.7
    axs[0].bar(range(len(circuits)), kl_tc, width=width)
    axs[0].set_title('TK_KL')
    for i,v in enumerate(kl_tc):
        axs[0].text(i, v + 1, str(v), color='black', fontweight='bold', fontsize='small',ha='center')
    axs[0].set_xticks(range(len(circuits)))
    axs[0].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    axs[1].bar(range(len(circuits)), sp_tc, width=width)
    axs[1].set_title('TK_SP')
    for i,v in enumerate(sp_tc):
        axs[1].text(i, v + 1, str(v), color='black', fontweight='bold', fontsize='small',ha='center')
    axs[1].set_xticks(range(len(circuits)))
    axs[1].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    axs[2].bar(range(len(circuits)), ga_kl_tc, width=width)
    axs[2].set_title('TK_GA_KL')
    for i,v in enumerate(ga_kl_tc):
        axs[2].text(i, v + 1, str(v), color='black', fontweight='bold', fontsize='small',ha='center')
    axs[2].set_xticks(range(len(circuits)))
    axs[2].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    axs[3].bar(range(len(circuits)), ga_sp_tc, width=width)
    axs[3].set_title('TK_GA_SP')
    for i,v in enumerate(ga_sp_tc):
        axs[3].text(i, v + 1, str(v), color='black', fontweight='bold', fontsize='small',ha='center')
    axs[3].set_xticks(range(len(circuits)))
    axs[3].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    axs[4].bar(range(len(circuits)),ga_tc,width=width)
    axs[4].set_title('TK_GA')
    for i,v in enumerate(ga_tc):
        axs[4].text(i,v + 1,str(v),color ='black', fontweight='bold', fontsize='small',ha='center')
    axs[4].set_xticks(range(len(circuits)))
    axs[4].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    plt.show()

def visualize_tc_improvement_to_Zomorodi():
    circuits = ['Sym9_147','4gt5_76','4mod7','Figure4','Parity_247','QFT4','QFT8','QFT16']
    kl_tc =    [70.84,57.14,60,50,0,50,68.42,75.94]
    sp_tc =    [62.5,28.57,60,50,0,50,78.95,80.45]
    ga_tc =    [75,71.43,60,50,0,50,78.95,87.97]
    ga_kl_tc = [87.5,71.43,60,50,0,50,78.95,54.89]
    ga_sp_tc = [75,71.43,60,50,0,50,78.95,54.89]

    fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(30,10))
    fig.suptitle('Verbesserung der Teleportationskosten in Prozent')

    width = 0.5
    axs[0].bar(range(len(circuits)), kl_tc, width=width)
    axs[0].set_title('TK_KL')
    for i,v in enumerate(kl_tc):
        axs[0].text(i, v + 1, str(v), color='black', fontweight='bold', fontsize='small',ha='center')
    axs[0].set_xticks(range(len(circuits)))
    axs[0].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    axs[1].bar(range(len(circuits)), sp_tc, width=width)
    axs[1].set_title('TK_SP')
    for i,v in enumerate(sp_tc):
        axs[1].text(i, v + 1, str(v), color='black', fontweight='bold', fontsize='small',ha='center')
    axs[1].set_xticks(range(len(circuits)))
    axs[1].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    axs[2].bar(range(len(circuits)), ga_kl_tc, width=width)
    axs[2].set_title('TK_GA_KL')
    for i,v in enumerate(ga_kl_tc):
        axs[2].text(i, v + 1, str(v), color='black', fontweight='bold', fontsize='small',ha='center')
    axs[2].set_xticks(range(len(circuits)))
    axs[2].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    axs[3].bar(range(len(circuits)), ga_sp_tc, width=width)
    axs[3].set_title('TK_GA_SP')
    for i,v in enumerate(ga_sp_tc):
        axs[3].text(i, v + 1, str(v), color='black', fontweight='bold', fontsize='small',ha='center')
    axs[3].set_xticks(range(len(circuits)))
    axs[3].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    axs[4].bar(range(len(circuits)),ga_tc,width=width)
    axs[4].set_title('TK_GA')
    for i,v in enumerate(ga_tc):
        axs[4].text(i,v + 1,str(v),color ='black', fontweight='bold', fontsize='small',ha='center')
    axs[4].set_xticks(range(len(circuits)))
    axs[4].set_xticklabels(circuits, rotation=65, ha='right', fontsize='small')

    plt.show()


#visualize_teleportation_cost_k_equals_2()
#visualize_teleportation_cost_k_equals_3()
#visualize_teleportation_cost_k_equals_4()
visualize_tc_improvement_to_Zomorodi()
