from ipywidgets import interact, interact_manual, widgets, Label
import numpy as np
import functools
from time import time, sleep
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

  lambda_T4=34
  v_sum=np.zeros((len_v_bits-lambda_T4,),dtype=np.uint64)
  for i in range(lambda_T4):
    v_sum+=v_bits[i:len_v_bits-lambda_T4+i]

  T4_PS=((v_sum==0) | (v_sum==lambda_T4)).sum()

  v_T5=np.zeros((5000,),dtype=np.uint64)
  for tau in range(5000):
    v_T5[tau]=(v_bits[:5000]^v_bits[tau+1:tau+5000+1]).sum()
  l_to_ret=(T1,T2,[[m_T3[0,0],m_T3[0,1]],[m_T3[1,0],m_T3[1,1]],[m_T3[2,0],m_T3[2,1]],[m_T3[3,0],m_T3[3,1]],[m_T3[4,0],m_T3[4,1]],[m_T3[5:,0].sum(),m_T3[5:,1].sum()]],T4_PS,v_T5)
  print(l_to_ret)
  return l_to_ret

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
  
  
raw_bytes=widgets.Textarea(value='',placeholder='',description='',disabled=True,layout=widgets.Layout(width='200px',height='200px'))
raw_bytes_label=widgets.HTML(value='<p><b>raw file</b></p>',placeholder='',description='',disabled=True, continuous_update=True,layout=widgets.Layout(width='150px',height='50px'))
comp_bytes=widgets.Textarea(value='',placeholder='',description='',disabled=True,layout=widgets.Layout(width='200px',height='200px'))
comp_bytes_label=widgets.HTML(value='<p><b>compressed file</b></p>',placeholder='',description='',disabled=True, continuous_update=True,layout=widgets.Layout(width='150px',height='50px'))
tests_results=widgets.HTML(value='',placeholder='',description='',disabled=True, continuous_update=True,layout=widgets.Layout(width='200px',height='300px'))

#entropy_value_widgets=widgets.FloatSlider(value=0.5,min=0,max=1,step=0.01,description='',disabled=False, continuous_update=False,orientation='horizontal',readout=True, readout_format='.1f')

nb_p=10000
v_p=0.5*np.arange(1,nb_p+1,dtype=np.float64)/(nb_p+1)
v_h=-1*(np.log2(v_p)*(v_p)+np.log2(1-v_p)*(1-v_p))

def gzipVsentropy(entropy_in=0.5):
  s_color_code_succeed='00FF00'
  s_color_code_failed='FF0000'

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

  v_in=generate_biased_rnd_vector_bits(p_target,2500)
  v_comp=comp_gzip(v_in)
  v_comp_bin=np.zeros((v_comp.shape[0]*8,),dtype=np.uint8)
  for i in range(8):
    v_comp_bin[i::8]=(v_comp>>i)&0x01

  s_line_raw=''
  s_line_comp=''  
  for i in range(v_in.shape[0]):
    s_line_raw='%s%d'%(s_line_raw,v_in[i])
  for i in range(v_comp_bin.shape[0]):
    s_line_comp='%s%d'%(s_line_comp,v_comp_bin[i])

  
  raw_bytes.value=s_line_raw
  comp_bytes.value=s_line_comp
  (T1,T2,M_T3,T4_PS,v_T5)=run_fips_tests(v_in)
  s_label_test='<table><tr><td>raw file size </td><td>&nbsp</td><td>%d bytes</td></tr><tr><td>compressed file size </td><td>&nbsp</td><td>%d bytes</td></tr></table>'%(int(v_in.shape[0]/8),int(v_comp_bin.shape[0]/8))
  
  if abs(T1-v_in.shape[0]/2)<346:
    s_color_code=s_color_code_succeed
  else:
    s_color_code=s_color_code_failed
  s_label_test='%s<table><tr><td style="color:#%s";>T1=%d</td><td>&nbsp&nbsp</td>'%(s_label_test,s_color_code,T1)

  if (T2<57.4)&(T2>1.03):
    s_color_code=s_color_code_succeed
  else:
    s_color_code=s_color_code_failed
  s_label_test='%s<td style="color:#%s";>T2=%.1f</td><td>&nbsp&nbsp</td>'%(s_label_test,s_color_code,T2)

  m_T3_bounds=[[2267,2733],[1079,1421],[502,748],[223,402],[90,223],[90,223]]

  s_label_T3='<table>'
  s_label_T3='%s<tr>'%(s_label_T3)
  s_label_T3="%s <td>&nbsp&nbsp\u03BB</td><td>&nbsp&nbsp</td>"%(s_label_T3)
  for i in range(6):
    s_label_T3="%s <td>&nbsp&nbsp%d</td><td>&nbsp&nbsp</td>"%(s_label_T3,i)
  s_label_T3='%s</tr>'%(s_label_T3)

  for b in range(2):
    s_label_T3='%s<tr>'%(s_label_T3)
    s_label_T3="%s <td>T<sub>3</sub>(%d,\u03BB)</td><td>&nbsp&nbsp</td>"%(s_label_T3,b)
    for i in range(6):
      if (M_T3[i][b]<m_T3_bounds[i][1])&(M_T3[i][b]>m_T3_bounds[i][0]):
        s_color_code=s_color_code_succeed
      else:
        s_color_code=s_color_code_failed       
      s_label_T3='%s <td style="color:#%s";>%d</td><td>&nbsp&nbsp</td>'%(s_label_T3,s_color_code,M_T3[i][b])
    s_label_T3='%s</tr>'%(s_label_T3)
  s_label_T3='%s</table>'%(s_label_T3)    

  if T4_PS==False:
    s_color_code=s_color_code_succeed
    s_result='SUCCEED'
  else:
    s_color_code=s_color_code_failed
    s_result='FAIL'
  s_label_test='%s<td style="color:#%s";>T4 %s</td><td>&nbsp&nbsp</td>'%(s_label_test,s_color_code,s_result)

  T5=v_T5.max()
  if (T5<2674)&(T5>2326):
    s_color_code=s_color_code_succeed
  else:
    s_color_code=s_color_code_failed

  s_label_test='%s<td style="color:#%s";>T5=%d</td></tr></table>'%(s_label_test,s_color_code,T5)

  s_label_test='%s<p></p>%s'%(s_label_test,s_label_T3)

  tests_results.value=s_label_test
  

#entropy_value_widgets.on_trait_change(generate_rnd, remove=False)

im=interact(gzipVsentropy)
im.widget.close()
im.widget.children[0].max=1.0
im.widget.children[0].min=0.0
im.widget.children[0].step=0.001
im.widget.children[0].description='Entropy'
im.widget.children[0].readout_format='.4f'
v_box_top=widgets.VBox([im.widget.children[0],widgets.HBox([widgets.VBox([raw_bytes,raw_bytes_label]),widgets.VBox([comp_bytes,comp_bytes_label]),tests_results])])
