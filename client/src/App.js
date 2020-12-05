import React from "react";
import { BrowserRouter as Router, Switch, Route } from "react-router-dom";
import SearchPage from "./SearchPage"
function App() {
  return (
    <Router>
      <Switch>
        <Route path="/" component={SearchPage} />
        <Route path="*" component={SearchPage} />
        {/* <Route
          exact
          path="/results"
          component={(props) => <Results {...props} />}
        /> */}
      </Switch>
    </Router>
  );
}

export default App;
