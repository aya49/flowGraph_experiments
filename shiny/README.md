the website is created using shiny, and styled using the codyhouse framework.

to install the codyhouse framework:

```
# first install npm

sudo apt update
sudo apt install nodejs
nodejs --version # verify

# install development tools

sudo apt install build-essential
# sudo apt remove nodejs npm # if you want to uninstall Node.js and npm packages

## Gulp configuration

# install the modules the framework requires for compiling SCSS into CSS
npm install

# launch your project on a development server
npm run gulp watch
```