#!/usr/bin/env python
# coding: utf-8

# # New Data Files: $\pi^+ , \pi^-, \pi^0 $

# Creating the separate dataframes.

# In[1]:


import pandas as pd

df_plus = pd.read_csv("pi_plus.csv")


# In[27]:


df_plus


# In[3]:


df_minus = pd.read_csv("pi_minus.csv")


# In[4]:


df_minus


# In[5]:


df_zero = pd.read_csv("pi_zero.csv")


# In[6]:


df_zero


# ## Splitting a dataframe.

# In[28]:


df_plus_quarter = df_plus.sample(frac=0.5)
df_plus_quarter


# In[29]:


df_minus_quarter = df_minus.sample(frac=0.5)
df_minus_quarter


# In[30]:


df_zero_half = df_zero.sample(frac = 1)
df_zero_half


# ## Combining multiple dataframes.

# In[31]:


df_combo = [df_plus_quarter, df_minus_quarter, df_zero_half]
together = pd.concat(df_combo)
together


# ## Creating Testing and Training Data

# In[32]:


train = together.sample(frac=0.75)
test = together.sample(frac=0.25)

train.to_csv("train.csv")
test.to_csv("test.csv")

