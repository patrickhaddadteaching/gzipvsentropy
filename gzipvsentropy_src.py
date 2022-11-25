from ipywidgets import widgets,interactive
import numpy as np
import gzip


def run_fips_tests(v_bits):
    len_v_bits=v_bits.shape[0]
    nb_nibble=int((v_bits.shape[0])/4)
    v_nibble=np.zeros((nb_nibble,),dtype=np.uint8)
    for i in range(4):
        v_nibble+=np.array(v_bits[i:4*nb_nibble:4],dtype=np.uint8)*(2**i)

    T1=v_bits.sum()

    T2=0
    for j in range(16):
        fj=(v_nibble==j).sum()
        T2+=((fj-5000.0/16)**2)/(5000.0/16)

    nb_lambda=32
    m_T3=np.zeros((nb_lambda,2),dtype=np.float64)
    for lambda_i in range(nb_lambda):
        word_size_in_bits=lambda_i+3

        nb_words=int(len_v_bits-word_size_in_bits)
        v_words=np.zeros((nb_words,),dtype=np.uint64)
        for i in range(word_size_in_bits):
            v_words+=np.array(v_bits[i:i+nb_words],dtype=np.uint8)*(2**i)    

        patern_0=2**(word_size_in_bits-1)+1
        patern_1=(2**(word_size_in_bits)-2**(word_size_in_bits-1)-2)
        m_T3[lambda_i,0]=(v_words==patern_0).sum()
        m_T3[lambda_i,1]=(v_words==patern_1).sum()

    v_cnt=np.zeros((2,),dtype=np.uint64)
    cnt_tmp=0
    for i in range(len_v_bits-1):
        if v_bits[i]==v_bits[i+1]:
            cnt_tmp+=1
            if v_cnt[int(v_bits[i+1])]<cnt_tmp:
                v_cnt[int(v_bits[i+1])]=cnt_tmp
        else:
            cnt_tmp=1
                
    v_T4=v_cnt
    
    v_T5=np.zeros((5000,),dtype=np.uint64)
    for tau in range(5000):
        v_T5[tau]=(v_bits[:5000]^v_bits[tau+1:tau+5000+1]).sum()
    return (T1,T2,[[m_T3[0,0],m_T3[0,1]],[m_T3[1,0],m_T3[1,1]],[m_T3[2,0],m_T3[2,1]],[m_T3[3,0],m_T3[3,1]],[m_T3[4,0],m_T3[4,1]],[m_T3[5:,0].sum(),m_T3[5:,1].sum()]],v_T4,v_T5)

image_height=155
image_width=150
param_height=15
cell_height=30
cell_width=100
cell_width_col0=130
cell_width_col1=100
cell_width_collast=150
cell_width_title=cell_width_col0+cell_width_col1+cell_width_collast+cell_width*2

nb_p=10000
v_p=0.5*np.arange(1,nb_p+1,dtype=np.float64)/(nb_p+1)
v_h=-1*(np.log2(v_p)*(v_p)+np.log2(1-v_p)*(1-v_p))


box_layout = widgets.Layout(height='%dpx'%cell_height,width='%dpx'%cell_width,border='1px solid gray',margin='0px 0px 0px 0px',justify_content='center',align_items='center')
box_layout2 = widgets.Layout(height='%dpx'%cell_height,width='%dpx'%(2*cell_width),border='1px solid gray',margin='0px 0px 0px 0px',justify_content='center',align_items='center')
box_layout_col0 = widgets.Layout(height='%dpx'%cell_height,width='%dpx'%(cell_width_col0),border='1px solid gray',margin='0px 0px 0px 0px',justify_content='center',align_items='center')
box_layout_col1 = widgets.Layout(height='%dpx'%cell_height,width='%dpx'%(cell_width_col1),border='1px solid gray',margin='0px 0px 0px 0px',justify_content='center',align_items='center')
box_layout_collast = widgets.Layout(height='%dpx'%cell_height,width='%dpx'%(cell_width_collast),border='1px solid gray',margin='0px 0px 0px 0px',justify_content='center',align_items='center')
box_layout_title = widgets.Layout(height='%dpx'%cell_height,width='%dpx'%(cell_width_title),border='1px solid gray',margin='0px 0px 0px 0px',justify_content='center',align_items='center')
box_layout_down = widgets.Layout(height='%dpx'%cell_height,width='%dpx'%(cell_width_title),margin='0px 0px 0px 0px',justify_content='center',align_items='center')

v_items=[]

v_items_h=[]
v_items_h.append(widgets.HTML(value='<p align="center"><b>Procedure A of AIS31</b></p>',layout=box_layout_title))
v_items.append(v_items_h)

v_items_h=[]
v_items_h.append(widgets.HTML(value='<p align="center">Statistical Test</p>',layout=box_layout_col0))
v_items_h.append(widgets.HTML('<p align="center">Test Result</p>',layout=box_layout_col1))
v_items_h.append(widgets.HTML('<p align="center">Test value</p>',layout=box_layout2))
v_items_h.append(widgets.HTML('<p align="center">Valid range</p>',layout=box_layout_collast))
v_items.append(v_items_h)

v_items_h=[]
v_items_h.append(widgets.HTML('<p align="center">Monobit test</p>',layout=box_layout_col0))
v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout_col1))
v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout2))
v_items_h.append(widgets.HTML('<p align="center">[9655 ; 10345]</p>',layout=box_layout_collast))
v_items.append(v_items_h)

v_items_h=[]
v_items_h.append(widgets.HTML('<p align="center">Poker test</p>',layout=box_layout_col0))
v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout_col1))
v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout2))
v_items_h.append(widgets.HTML('<p align="center">[1.03 ; 57.4]</p>',layout=box_layout_collast))
v_items.append(v_items_h)

m_T3_bounds=[[2267,2733],[1079,1421],[502,748],[223,402],[90,223],[90,223]]
for i in range(6):
    v_items_h=[]
    if i==5:
        v_items_h.append(widgets.HTML('<p align="center">Run test &#x3BB &ge; %d</p>'%(i+1),layout=box_layout_col0))
    else:
        v_items_h.append(widgets.HTML('<p align="center">Run test &#x3BB=%d</p>'%(i+1),layout=box_layout_col0))
    v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout_col1))
    v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout))
    v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout))
    v_items_h.append(widgets.HTML('<p align="center">[%d ; %d]</p>'%(m_T3_bounds[i][0],m_T3_bounds[i][1]),layout=box_layout_collast))
    v_items.append(v_items_h)

v_items_h=[]
v_items_h.append(widgets.HTML('<p align="center">Long Run test</p>',layout=box_layout_col0))
v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout_col1))
v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout2))
v_items_h.append(widgets.HTML('<p align="center">[0 ; 33]</p>',layout=box_layout_collast))
v_items.append(v_items_h)


v_items_h=[]
v_items_h.append(widgets.HTML('<p align="center">Autocorrelations test</p>',layout=box_layout_col0))
v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout_col1))
v_items_h.append(widgets.HTML('<p align="center"> </p>',layout=box_layout2))
v_items_h.append(widgets.HTML('<p align="center">[2326 ; 2674] </p>',layout=box_layout_collast))
v_items.append(v_items_h)

nb_lines=len(v_items)
v_hb=[]
for i in range(nb_lines):
    v_tmp=v_items[i]
    nb_row=len(v_tmp)
    v_box=[]
    for j in range(nb_row):
        v_box.append(v_tmp[j])
    v_hb.append(widgets.HBox(v_box))
v_box_table=widgets.VBox(v_hb)


entropy_value_widgets=widgets.FloatSlider(description='RAW BITS ENTROPY',style = {'description_width': 'initial'},min=0,max=1.0,step=0.001,layout=widgets.Layout(height='%dpx'%(2*param_height),width='%dpx'%(2*image_width+cell_width_title)))
raw_bytes=widgets.Textarea(value='',rows=5000,placeholder='',description='',disabled=True,layout=widgets.Layout(height='%dpx'%(image_height),width='%dpx'%(2*image_width)))
raw_bytes_comp=widgets.Textarea(value='',rows=5000,placeholder='',description='',disabled=True,layout=widgets.Layout(height='%dpx'%(image_height),width='%dpx'%(2*image_width)))


def comp_gzip(v_bits_in):
    n_in=v_bits_in.shape[0]
    v_byte_to_int=np.zeros((int(n_in/8),),dtype=np.uint8)
    for i in range(8):
        v_byte_to_int+=(v_bits_in[i::8]<<i)  
    v_out = np.frombuffer(gzip.compress(v_byte_to_int),dtype=np.uint8)
    return v_out

def generate_biased_rnd_vector_bits(p_in,n_in):
    v_byte_to_ret=np.zeros((n_in,),dtype=np.uint8)
    v_raw_rnd_bin=np.array(np.random.rand(8*n_in)<p_in,dtype=np.uint8)
    return v_raw_rnd_bin

def generate_rnd(entropy_in):
    
    

    if entropy_in==1:
        p_target=0.5
    elif entropy_in==0:
        p_target=np.random.randint(0,2)
    else:
        rnd=np.random.randint(0,2)
        if rnd==0:
            p_target=v_p[np.where(v_h>entropy_in)[0][0]]
        else:
            p_target=1-v_p[np.where(v_h>entropy_in)[0][0]]   
    
    v_raw_bits=generate_biased_rnd_vector_bits(p_target,2500)
    nb_bits=v_raw_bits.shape[0]
    v_comp=comp_gzip(v_raw_bits)
    v_comp_bin=np.zeros((v_comp.shape[0]*8,),dtype=np.uint8)
    for i in range(8):
        v_comp_bin[i::8]=(v_comp>>i)&0x01 
        
    nb_nibble=int(v_raw_bits.shape[0]/4)
    v_nibble=np.zeros((nb_nibble,),dtype=np.uint8)
    for i in range(4):
        v_nibble+=np.array(v_raw_bits[i:4*nb_nibble:4],dtype=np.uint8)*(2**i)
    s_line_raw_bytes_value='RAW DATA : %d BITS\n'%(nb_bits)
    for i in range(nb_nibble):
        s_line_raw_bytes_value='%s%X'%(s_line_raw_bytes_value,v_nibble[i])
    raw_bytes.value=s_line_raw_bytes_value 

    nb_nibble_comp=int(v_comp_bin.shape[0]/4)
    v_nibble_comp=np.zeros((nb_nibble_comp,),dtype=np.uint8)
    for i in range(4):
        v_nibble_comp+=np.array(v_comp_bin[i:4*nb_nibble_comp:4],dtype=np.uint8)*(2**i)
    s_line_comp_bytes_value=''
    s_line_comp_bytes_value='COMPRESSED DATA : %d BITS\n'%(v_comp_bin.shape[0])
    for i in range(nb_nibble_comp):
        s_line_comp_bytes_value='%s%X'%(s_line_comp_bytes_value,v_nibble_comp[i])
    raw_bytes_comp.value=s_line_comp_bytes_value 
    

    (T1_o,T2_o,M_T3_o,v_T4_o,v_T5_o)=run_fips_tests(v_raw_bits)
    
    v_items[2][2].value='<p align="center">%d</p>'%T1_o
    if abs(T1_o-nb_bits/2)<346:
        v_items[2][1].value='<p align="center"; style="color:#00FF00";> PASS </p>'
    else:
        v_items[2][1].value='<p align="center"; style="color:#FF0000";> FAIL </p>'    
    
    v_items[3][2].value='<p align="center">%.2f</p>'%T2_o
    if (T2_o<57.4)&(T2_o>1.03):
        v_items[3][1].value='<p align="center"; style="color:#00FF00";> PASS </p>'
    else:
        v_items[3][1].value='<p align="center"; style="color:#FF0000";> FAIL </p>'    

    for i in range(6):
        cnt_error=0
        for b in range(2):
            v_items[4+i][2+b].value='<p align="center">%d</p>'%M_T3_o[i][b]
            if (M_T3_o[i][b]<m_T3_bounds[i][1])&(M_T3_o[i][b]>m_T3_bounds[i][0]):
                cnt_error+=1                
        if (cnt_error>0):
            v_items[4+i][1].value='<p align="center"; style="color:#00FF00";> PASS </p>'
        else:
            v_items[4+i][1].value='<p align="center"; style="color:#FF0000";> FAIL </p>'
    
    longest_run=v_T4_o.max()
    v_items[10][2].value='<p align="center">Longest run : %d</p>'%(longest_run)
    if (longest_run<34):
        v_items[10][1].value='<p align="center"; style="color:#00FF00";> PASS </p>'
    else:
        v_items[10][1].value='<p align="center"; style="color:#FF0000";> FAIL </p>'   

    max_autocor=v_T5_o.max()
    min_autocor=v_T5_o.min()
    v_items[11][2].value='<p align="center">MIN=%d MAX=%d</p>'%(min_autocor,max_autocor)
    if (max_autocor<2674)&(max_autocor>2326)&(min_autocor<2674)&(min_autocor>2326):
        v_items[11][1].value='<p align="center"; style="color:#00FF00";> PASS </p>'
    else:
        v_items[11][1].value='<p align="center"; style="color:#FF0000";> FAIL </p>'  
        
interactive_plot = interactive(generate_rnd,{'manual': True}, entropy_in =entropy_value_widgets)
interactive_plot.children[-2].description='GENERATE RANDOM NUMBERS'
interactive_plot.children[-2].layout=widgets.Layout(width='%dpx'%(2*image_width))  
v_box_top=widgets.VBox([interactive_plot.children[0],widgets.HBox([widgets.HBox([widgets.VBox([interactive_plot.children[-2],raw_bytes,raw_bytes_comp])]),v_box_table])])
        
