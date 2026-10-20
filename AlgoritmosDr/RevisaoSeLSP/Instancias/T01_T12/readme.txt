practical problems:	12 problem files with 8 periods and 2 machines

The structure of the 12 problem files (p_8p2l.??0) is as follows:
-------------------------------------------------
maximal number of products / number of machines ("# Anzahl Produkte", not all of the products are used in the current problem)
number of macro periods ("# Anzahl Makroperioden")
maximal number of micro periods per macro period ("# Anzahl Mikroper. je Makroper.", not all of them need to be used)

location of files for the following data:
-----------------------------------------
initial setup state per machine ("..\data\8per2lin\y_8p2l.??0", 0 is a neutral initial state) 
production costs per machine and product ("...c_...")
changeover costs per machine for changeover from product i to j ("...s_...", "#Zusatzprodukt" for changeover from neutral state 0 to another product per machine)
inventory holding costs per unit of product and macro period ("...h_...")
demand per product and macro period ("...d_...")
minimum lotsize per product("...m_...")
capacity consumption per line and unit of product ("...a_...") 
capacity per line and macro period ("...k_...")
initial inventory per product ("...i_...")
changeover times per machine for changeover from product i to j ("...st_...", "#Zusatzprodukt" for changeover from neutral state 0 to another product per machine)

