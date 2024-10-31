# VoxSoft

Welcome to the **VoxSoft**, a Unity project that enables anyone to create volumetric meshes of soft robots using simple shape functions, simulate them using Extended Position-Based Dynamics (XPBD), and optimize their designs for better performance.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Usage](#usage)
  - [Creating Volumetric Meshes](#creating-volumetric-meshes)
  - [Simulating with XPBD](#simulating-with-xpbd)
  - [Design Optimization](#design-optimization)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Introduction

Soft robotics is a rapidly growing field focusing on creating robots with flexible and compliant materials. This project provides a user-friendly platform to design, simulate, and optimize soft robots within Unity.

## Features

- **Volumetric Mesh Generation**: Create complex soft robot geometries using simple shape functions like spheres, cylinders, and custom shapes.
- **Extended Position-Based Dynamics Simulation**: Utilize XPBD for realistic and stable simulation of soft body physics.
- **Design Optimization**: Optimize robot designs based on performance metrics such as movement efficiency, material usage, or response time.
- **User-Friendly Interface**: Intuitive tools and interfaces suitable for both beginners and advanced users.
- **Extensibility**: Modular codebase allowing for easy customization and integration with other systems.

## Getting Started

### Prerequisites

- **Unity Engine**: Version 2022.3.13f1 or later.
- **Git**: For cloning the repository (optional if downloading ZIP).

### Installation

1. **Clone the Repository**

   ```bash
   git clone https://github.com/yourusername/soft-robot-designer.git

2.	Open the Project in Unity
	•	Launch Unity Hub.
	•	Click “Open” and select the cloned project folder.
	•	Allow Unity to import all assets.

Usage

Open the VoxSoftSim Scene and enable the Soft Body Controller. Default: 
Num Sub Steps = 50
Edge Compliance = 0.1
Vol Compliance = 0
Damping Coefficient = 0.99
Pressure = 0
Scale = 0.005


Acknowledgments

	•	Unity Technologies for the game engine.
	•	Researchers and developers who contributed to the development of XPBD.
	•	The soft robotics community for inspiration and resources.
	•	Habrador for his initial implementation of XPBD in Unity (https://github.com/Habrador/Ten-Minute-Physics-Unity)

Feel free to open an issue or submit a pull request for suggestions and improvements!
