# LCP Simulator
This repository contains all underlying scripts of the LCP Simulator and a description of their implementation in the web application builder Anvil Works. Details on what the LCP Simulator does and how to use this web tool can be found in detail on the website itself which is equipped with extensive guidelines and help and error messages.

The Linear CoPolymer Simulator (LCP Simulator) is a Python-based, user-friendly web tool available to anyone at: https://lcp-simulator.anvil.app. It can simulate linear polymers consisting of two different building blocks (A and B) and their behavior during (enzymatic) cleavage and/or analytics, and is supposed to be of use for anyone performing analytics or molecular biology on those type of copolymers. 

# Implementation of Python code in Anvil Works
To provide our set of Python scripts for simulating linear copolymer analytics and cleavage with a
comprehensible graphical user interface (GUI), and to make it easily accessible as a web tool, we used
the commercial web application builder Anvil Works for implementation. Anvil Works hosts the
website and offers the opportunity to design Python-based GUIs from predefined building blocks.
There is a distinction between two types of code: First, there is server code, which includes the
computationally intensive functions that form the true core of what the simulations are designed to
compute. Second, there is client code that defines the behavior of the GUI, triggers events based on
user clicks, and calls functions from the server code.

## Server code
In our case, the server code consists of just one module called server_module.py within the
supplementary data. Basically, it is a collection of two types of functions, simple functions and callable
functions. Callable functions have the decorator “@anvil.server.callable” and can thus be called from
the client code. Simple functions do not have this decorator, so they are not callable from the client
code, but can be called from all other functions within the server code. Whereas you cannot use
Python packages like pandas, numpy, or matplotlib within the client code, you can use them within
the server code.

## Client code
Our client code is composed of six forms, one for each of the five tabs and an additional one for the
overall framework. The framework contains the sidebar of the website with its five tabs and controls
which of the other five forms is visible and accessible depending on which tab the user selects.
In contrast to server code modules which are code only, client code forms are split into a design tab
and a code tab. In the former, features such as text boxes, buttons, or cards are added by drag and
drop, which can then be modified, for example, in terms of color and style, or by adding tooltips. In
the latter, Python code is entered which controls what happens when design elements are clicked or
filled in by the user of the web tool. For example, a click on the “Submit” button created in the design
tab leads to a call of the corresponding computationally intensive core function from the server code,
as defined within the code tab, which ultimately leads to the results that form the output for the user.
The supplementary data includes the contents of the code tab for each of the six client code forms,
those for the five tabs (home.py, DP_histogram.py, product_profile.py, pattern_analysis.py,
block_sizes.py) and the one for the framework (app_frame.py).
